/*    Copyright@Louxin Zhang  National University of Singapore
*
* *  This program takes a set of rooted binary phylogenetic 
*     trees and compute tree-child networks (acyclic rooted networks with
*     a set of leaves satisfying some constraints)
* *   that contain every input tree and have the minimum  number of reticulations,
*      i.e., no networks with less reticulations can ontain all the trees.
* *
*    Input tree files contain a set of trees given in Newick format.
* *  We also assume that trees nodes are represented by integers. 
*    For example, the following is a file containing 5 trees with 20 taxa each
*    ------------  begin of the file ----------------
*  (((((10,14),20),1),((18,7),(15,(13,9)))),((2,(5,12)),((((8,17),(16,3)),11),((4,19),6))))
*  (((2,(5,12)),15),((((10,14),20),(13,9)),((11,((16,1),3)),(((4,19),((18,(8,17)),7)),6))))
*  ((6,(18,(13,9))),(((2,(5,12)),15),((1,(((8,17),((16,3),7)),11)),(((10,14),20),(4,19)))))
*  ((((6,(2,(5,12))),((10,14),20)),((18,(8,17)),(15,(13,9)))),(((((16,1),3),7),11),(4,19)))
*  (((11,1),(((18,(8,17)),7),(((10,14),20),(13,9)))),((6,((2,(5,12)),15)),((16,3),(4,19))))
*    ------------- end of the file -----------------
* *
* *   Compiling command:  gcc Arb_Ntk_Approx_NW.c -o AAST
* *   Run command: ./AAST <tree_input_file> <no_taxa> <output_file>
* *
*     The program will output networks into <output_file> if there are solutions. 
*     Limitations:  the no_taxa can be  as large as 128;
*     Two program constants CUT_OFF_SIZE and MEMORY_LIMIT control the runing time 
*
* */


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>

#define TREE_OR_LEAF 0
#define RETNODE 1
#define T_EDGE 'a'
#define RED_EDGE 'b'
#define T_LEAF -1
#define T_ROOT 2
#define T_OTHER 1

#define MAX_DEG 3
#define NO_TREES 100
#define NO_NODES  260
#define NO_LEAVES 128
#define MAXSIZE  1024   /* deal with ntks with at most 20 leaves */
#define TWO 2
#define  THREE 3
#define TEST 4   /* test purpose */
#define DIGIT_NO 100 

#define CUT_OFF_SIZE 4000
#define MEMORY_LIMIT 500000;



/* encode an edge (x, y) by a short integer: z*100+y to save memory space */
typedef struct graph {
    short *left_ends;    
    short *right_ends;    
    short  *parent[NO_NODES];    /* two parents p1*100+p3 */
    short no_nodes;
    short no_edges;
    struct graph *next;
} g_set;

/* a graph is represented by a linked list of edges (x, y), 
 * which are encoded by an interger using * ENCODE(1000)*x+y */
typedef struct f_gph {
   /* short encoded_edge[MAXSIZE]; */
   int *lt_ends;
   int *rt_ends;
   short no_edges;
   struct f_gph *next;
} ntks_struct;




struct node {
	char *item;
        short  score;
	struct node* next;
};

struct node *Q0_front, *Q1_front;
struct node *Q0_end, *Q1_end;


void Qinsert(struct node **front, struct node **end,  char *str) {
     struct node *ptr;

       /*  printf("Inside Qinsert\n"); */
	ptr = (struct node *) malloc (sizeof(struct node));
	if (ptr == NULL) {
		printf("\nOVERFLOW\n");
		return;
	} else {
		ptr ->item = (char *)malloc(strlen(str)+1);
                strcpy(ptr ->item, str); ptr->score=0; ptr ->next =NULL;
		if (*front == NULL) { *front = ptr; *end = ptr;
		} else { (*end) -> next = ptr; (*end) = ptr; }
	}
}

short  Qdelete(struct node **front, struct node **end, short *score,  char *str) {  
    struct node *ptr;


    if (*front == NULL) {  
        *score=0; return 0;  
    }  else {
        ptr=(*front); strcpy(str, ptr->item);  *score=ptr->score;
        if (*front==*end) { *front=NULL; *end=NULL;
        } else { (*front) = ptr -> next;  }
        free(ptr->item);  free(ptr);
        return 1;
    } 
}   
 
void Qremove(struct node **front, struct node **end){
   struct node *ptr;
   int k; 
    k=0;
    while (*front!=NULL) {
       k++; ptr=*front;
       if (*front==*end) {*front=NULL; *end=NULL;}
       else  (*front) = ptr -> next;
       free(ptr->item);  free(ptr);
    }
}

/* Merges two subarrays of arr[].  First subarray is arr[left..m]
  Second subarray is arr[m+1..right]
 */
 void Merge(short A[], short  left, short m, short right, short len[]) {
	short i, j, k;
        short  n1 = m - left + 1;
        short  n2 = right - m;
        short  L[n1], R[n2];

        for (i = 0; i < n1; i++) L[i] = A[left + i];
        for (j = 0; j < n2; j++) R[j] = A[m + 1 + j];

        i = 0; j=0; k = left; 
        while (i < n1 && j < n2) {
           if (len[L[i]] >= len[R[j]]) {
              A[k] = L[i]; i++;
           } else {
              A[k] = R[j]; j++;
           }
           k++;
       }

       while (i < n1) { A[k] = L[i]; i++; k++; }
       while (j < n2) { A[k] = R[j]; j++; k++; }
} /* Merge */


/* left is for left index and r is right index of the sub-array 
 of arr to be sorted 
*/
 void MergeSort(short  A[], short left, short right, short len[]) {
     short m;

    if (left < right) {
            m = left + (right - left) / 2;
            MergeSort(A, left, m, len);
            MergeSort(A, m + 1, right, len);
            Merge(A, left, m, right, len);
    }
}


void Count_NNN(struct node *A, struct node *B) {
  int i;
  struct node *ptr;

  ptr=A; i=0;
  while (ptr!=B) { i++; ptr=ptr->next; }
  i++;
  printf("size %d\n", i);
}

/* Sort the nodes of a ntk represented by ch[][ ] so that edges from a node to a node appearing later
 * in the list. It is easy to see that the first node of the list is the root and 
 * the last node must be a leaf */

void top_sort(short ch[][THREE], short no_nodes,
   short no_ntk_edges, short sorted_nodes[], short no_lfs){
    short i;
    short  *in_deg, *out_d;
    short  *roots;
    short  r_no, no, u;


    in_deg=(short *)malloc(sizeof(short)*(no_nodes+1));
    out_d=(short *)malloc(sizeof(short)*(no_nodes+1));
    roots=(short *)malloc(sizeof(short)*(no_nodes+1));

    in_deg[0]=0; out_d[0]=1;
    for (i=1; i<=no_lfs; i++) { in_deg[i]=1; out_d[i]=0; }

    no=2*no_lfs;
    for (i=no_lfs+2; i<=no_nodes; i++) {
       if (i<=no || i%2==1) { in_deg[i]=1; out_d[i]=2;}
       else { in_deg[i]=2; out_d[i]=1;}
    }


    roots[0]=0; r_no=1;
    no=0;
    while (r_no >0) {
      u=roots[r_no-1]; sorted_nodes[no]=u;
      no +=1; r_no=r_no-1;
      for (i=0; i<out_d[u]; i++) {
         in_deg[ch[u][i+1]]=in_deg[ch[u][i+1]]-1;
         if (in_deg[ch[u][i+1]]==0) {roots[r_no]=ch[u][i+1]; r_no +=1; }
      }
    } /* while loop */
    free(in_deg); free(out_d); free(roots);
}



/* convert a tree or a ntk that is kept as an array of edges into a child array and parent array */
/* Ch[i] keeps the children of the ith node. Ch[i][0] equals the number of childrens of
 *  the ith node. each tree node has two children, recorded in Ch[i][1] and Ch[i][2];
 *  the root and each ret have 1, kept in Ch[i][1].
 */
void Process_Tree(short tree[][TWO], short edges_no, short no_nodes,
 short ch[][THREE], short pa[]){
    short  i, k;

   for (i=0; i<=no_nodes;  i++) { ch[i][0]=0; pa[i]=-1;}
   for (i=1; i<=edges_no;  i++) {
      ch[tree[i][0]][0] +=1; k=ch[tree[i][0]][0]; ch[tree[i][0]][k]=tree[i][1];
      pa[tree[i][1]]=tree[i][0];
   }
} /* process100 */


void Compute_parents(short *gph[], short edges_no, short nodes_no, 
short pa[][NO_NODES]){
    short  i, m, n, k, min, pos;

   pa[0][0]=-1; 
   for (i=1; i<NO_NODES; i++) pa[i][0]=0;
   for (i=1; i<=edges_no;  i++) {
      k=pa[gph[1][i]][0];
      pa[gph[1][i]][0]=1+k;  pa[gph[1][i]][k+1]=gph[0][i];
   }
   /* sort the parents in the increasing order */
   for (i=1; i<=nodes_no; i++){
     k=pa[i][0];
     if (k>1) {
       for (m=1; m<=k; m++) {
          min=pa[i][m]; pos=m;
          for (n=m+1; n<=k; n++) {
             if (pa[i][n]< min) {
                min=pa[i][n]; pos=n;
             }
          }
          if (pos > m) { pa[i][pos]=pa[i][m]; pa[i][m]=min; }
       }
     }
   } /* for */

} /* Compute parent array */


short  check_edge(short a, short b,  short edges[][2], short no_edges){
      short i, j;

  for (i=0; i<no_edges; i++) { 
     if (a==edges[i][0] && b==edges[i][1]) { return 1;} 
  }
  return 0;
}



/* used */
/* gph is a tree */
short  Simplify_Ntks_Initialization(short gph[][2], short no_edges,
   short no_nodes, short leaf_map[], short current_comp,  short last){
    short i, j, k, r;
    short  in_deg[NO_NODES],ot_deg[NO_NODES]; 
    short  internal[NO_NODES],external[NO_NODES];
    short  v_type[NO_NODES];
    short no_lfs, no_trs_rets, cpy[2];

    for (i=0; i<NO_NODES; i++) { in_deg[i]=0; ot_deg[i]=0; }
    for (i=1; i<=no_edges; i++) { in_deg[gph[i][1]] +=1; ot_deg[gph[i][0]] +=1;}
    k=0; j=1;
    for  (i=0; i<NO_NODES; i++) {
       if (in_deg[i]==0 &&  ot_deg[i]==1) { v_type[i]=4;
        }  /* root */
       else if (in_deg[i]==1 &&  ot_deg[i]==2) {  
            v_type[i]=3; internal[i]=k; k +=1;  /* tree node */
       } else  if (in_deg[i]==1 &&  ot_deg[i]==1) {
                v_type[i]=2; /* degree 2*/
                internal[i]=k; k +=1;
       } else if (in_deg[i]==1 &&  ot_deg[i]==0) {
               v_type[i]=1;  /* leaf */
               external[i]=j; j +=1;
       } else if (in_deg[i]==2 &&  ot_deg[i]==1) {
                v_type[i]=0;  /* ret */
                internal[i]=k; k +=1;
       } else v_type[i]=-1;
    }

   no_lfs=j+1; no_trs_rets=k;
   for (i=1; i<=no_edges; i++) {
          if (v_type[gph[i][0]]==4) { r=i;   gph[i][0]=0; 
            leaf_map[0]=gph[i][1];
           } else if (v_type[gph[i][0]]==3 ||v_type[gph[i][0]]==0) 
              gph[i][0]= no_lfs + internal[gph[i][0]];

          if (v_type[gph[i][1]]==3 || v_type[gph[i][1]]==0)
              gph[i][1]= no_lfs + internal[gph[i][1]];
          else if (v_type[gph[i][1]]==1) { 
               leaf_map[external[gph[i][1]]]=gph[i][1];
               gph[i][1]=external[gph[i][1]]; 
          }
          /*
          printf("%d %d\n", gph[i][0],  gph[i][1]);
          */
  }

    cpy[0]=gph[1][0]; cpy[1]=gph[1][1];
    gph[1][0]=gph[r][0]; gph[1][1]=gph[r][1];
    gph[r][0]=cpy[0]; gph[r][1]=cpy[1];
    gph[0][0]=gph[no_edges][0]; gph[0][1]=gph[no_edges][1];
} /* simplify_process2 */

short  Simplify_Ntks(short gph[][2], short no_edges, short no_nodes,
    short leaf_map[]){
    short i, j, k, r;
    short  in_deg[NO_NODES],ot_deg[NO_NODES],internal[NO_NODES];
    short external[NO_NODES], v_type[NO_NODES], cpy[2];
    short no_lfs, no_trs_rets;

    for (i=0; i<NO_NODES; i++) { in_deg[i]=0; ot_deg[i]=0; }
    for (i=1; i<=no_edges; i++) { in_deg[gph[i][1]] +=1; ot_deg[gph[i][0]] +=1;}

    k=0; j=1;
    for  (i=0; i<NO_NODES; i++) {
       if (in_deg[i]==0 &&  ot_deg[i]==1) { v_type[i]=4; /* root */
       } else if (in_deg[i]==1 &&  ot_deg[i]==2) {  
            v_type[i]=3; internal[i]=k; k +=1;  /* tree node */
       } else  if (in_deg[i]==1 &&  ot_deg[i]==1) {
                v_type[i]=2; /* degree 2*/
                internal[i]=k; k +=1;
       } else if (in_deg[i]==1 &&  ot_deg[i]==0) {
               v_type[i]=1;  /* leaf */
               external[i]=j; j +=1;
       } else if (in_deg[i]==2 &&  ot_deg[i]==1) {
                v_type[i]=0;  /* ret */
                internal[i]=k; k +=1;
       } else v_type[i]=-1;
    }

   no_lfs=j+1; no_trs_rets=k; /* no_trs_rets is the number of internals nodes */
   for (i=1; i<=no_edges; i++) {
          if (v_type[gph[i][0]]==4) { 
              r=i;   gph[i][0]=0; 
              leaf_map[0]=no_lfs+internal[gph[i][1]];
           } else if (v_type[gph[i][0]]==3 ||v_type[gph[i][0]]==0) 
              /* tree nodes or ret */
              gph[i][0]= no_lfs + internal[gph[i][0]];

          if (v_type[gph[i][1]]==3 || v_type[gph[i][1]]==0)
              gph[i][1]= no_lfs + internal[gph[i][1]];
          else if (v_type[gph[i][1]]==1) { 
               /* gph[i][i] is a leaf, rename it as external[gph[i][1]] */ 
               leaf_map[external[gph[i][1]]]=gph[i][1];
               gph[i][1]=external[gph[i][1]]; 
          }
          /*
          printf("%d %d\n", gph[i][0],  gph[i][1]);
          */
    }

    cpy[0]=gph[1][0]; cpy[1]=gph[1][1]; /* copy the first edge  for swap */
    gph[1][0]=gph[r][0]; gph[1][1]=gph[r][1]; /* move the the edge leaving the root to the first edge */
    gph[r][0]=cpy[0]; gph[r][1]=cpy[1]; /* move the first edge to the r-th edge */
    gph[0][0]=gph[no_edges][0]; gph[0][1]=gph[no_edges][1]; /* copy the last edge to the 0-th edge */
} /* simplify_process2 */

short  Simplify_Ntks_July22(short *gph[], short no_edges){
    short i, j, k, r;
    short  *in_deg, *ot_deg;
    short *internal;
    short  *v_type, cpy[TWO];
    short no_lfs, no_trs_rets;
    int max_node_id;


    
    max_node_id=0;
    for (i=1; i<=no_edges; i++) { 
            if (gph[1][i]>max_node_id) max_node_id=gph[1][i];
            if (gph[0][i]>max_node_id) max_node_id=gph[0][i];
    }
   

    in_deg=(short*)malloc((max_node_id+2)*sizeof(short));
    ot_deg=(short*)malloc((max_node_id+2)*sizeof(short));
    internal=(short*)malloc((max_node_id+2)*sizeof(short));
    v_type=(short*)malloc((max_node_id+2)*sizeof(short));
    
    for (i=0; i<=max_node_id; i++) { in_deg[i]=0; ot_deg[i]=0; }

    for (i=1; i<=no_edges; i++) { 
          in_deg[gph[1][i]] +=1; ot_deg[gph[0][i]] +=1;
    }
    k=0; j=1;
    for  (i=0; i<=max_node_id; i++) {
       if (in_deg[i]==0 &&  ot_deg[i]==1) { v_type[i]=4;
        }  /* root */
       else if (in_deg[i]==1 &&  ot_deg[i]==2) {  
            v_type[i]=3; internal[i]=k; k +=1;  /* tree node */
       } else  if (in_deg[i]==1 &&  ot_deg[i]==1) {
                v_type[i]=2; /* degree 2*/
                internal[i]=k; k +=1;
       } else if (in_deg[i]==1 &&  ot_deg[i]==0) {
               v_type[i]=1;  /* leaf */
               j +=1; 
       } else if (in_deg[i]>=2 &&  ot_deg[i]==1) {
                v_type[i]=0;  /* ret */
                internal[i]=k; k +=1;
       } else v_type[i]=-1;
    }

   no_lfs=j+1; no_trs_rets=k;
   for (i=1; i<=no_edges; i++) {
          if (v_type[gph[0][i]]==4) { r=i;   gph[0][i]=0; 
              /*
              leaf_map[0]=gph[i][1];
              printf("------------hihi %d\n",  leaf_map[0]);
               */
           }
          if (v_type[gph[0][i]]==3 ||v_type[gph[0][i]]==0) 
              gph[0][i]= no_lfs + internal[gph[0][i]];
          if (v_type[gph[1][i]]==3 || v_type[gph[1][i]]==0)
              gph[1][i]= no_lfs + internal[gph[1][i]];

          if (v_type[gph[1][i]]==1) { 
               /*
               */
          }
    }
    cpy[0]=gph[0][1]; cpy[1]=gph[1][1];
    gph[0][1]=gph[0][r]; gph[1][1]=gph[1][r];
    gph[0][r]=cpy[0]; gph[1][r]=cpy[1];
    gph[0][0]=gph[0][no_edges]; gph[1][0]=gph[1][no_edges];

     /*
     printf("---- end 0: %d \n", no_edges); 
    for (i=1; i<=no_edges; i++) 
     printf("%d  %d\n", gph[0][i], gph[1][i]);
     printf("\n\n");
     */


    free(in_deg);
    free(ot_deg);
    free(internal);
    free(v_type);
} /* simplify_process2 */

short  Simplify_Ntks_Aug22(int gph[][2], short no_edges){
    int i;
    short  j, k, r;
    short  *in_deg, *ot_deg;
    short *internal;
    short  *v_type, cpy[TWO];
    short no_lfs, no_trs_rets;
    int max_node_id;

    max_node_id=0;
    for (i=1; i<=no_edges; i++) { 
            if (gph[i][1]>max_node_id) max_node_id=gph[i][1];
            if (gph[i][0]>max_node_id) max_node_id=gph[i][0];
    }
   
    /* printf("---- max: %d\n", max_node_id); */

    in_deg=(short*)malloc((max_node_id+1)*sizeof(short));
    ot_deg=(short*)malloc((max_node_id+1)*sizeof(short));
    internal=(short*)malloc((max_node_id+1)*sizeof(short));
    v_type=(short*)malloc((max_node_id+1)*sizeof(short));
    
    for (i=0; i<=max_node_id; i++) { in_deg[i]=0; ot_deg[i]=0; }

    for (i=1; i<=no_edges; i++) { in_deg[gph[i][1]] +=1; ot_deg[gph[i][0]] +=1;}
    k=0; j=1;
    for  (i=0; i<=max_node_id; i++) {
       if (in_deg[i]==0 &&  ot_deg[i]==1) { v_type[i]=4;
        }  /* root */
       else if (in_deg[i]==1 &&  ot_deg[i]==2) {  
            v_type[i]=3; internal[i]=k; k +=1;  /* tree node */
       } else  if (in_deg[i]==1 &&  ot_deg[i]==1) {
                v_type[i]=2; /* degree 2*/
                internal[i]=k; k +=1;
       } else if (in_deg[i]==1 &&  ot_deg[i]==0) {
               v_type[i]=1;  /* leaf */
               j +=1; 
       } else if (in_deg[i]>=2 &&  ot_deg[i]==1) {
                v_type[i]=0;  /* ret */
                internal[i]=k; k +=1;
       } else v_type[i]=-1;
    }

   no_lfs=j+1; no_trs_rets=k;
   for (i=1; i<=no_edges; i++) {
          if (v_type[gph[i][0]]==4) { r=i;   gph[i][0]=0; 
              /*
              leaf_map[0]=gph[i][1];
              printf("------------hihi %d\n",  leaf_map[0]);
               */
           }
          if (v_type[gph[i][0]]==3 ||v_type[gph[i][0]]==0) 
              gph[i][0]= no_lfs + internal[gph[i][0]];
          if (v_type[gph[i][1]]==3 || v_type[gph[i][1]]==0)
              gph[i][1]= no_lfs + internal[gph[i][1]];

          if (v_type[gph[i][1]]==1) { 
               /*
               leaf_map[external[gph[i][1]]]=gph[i][1];
               gph[i][1]=external[gph[i][1]]; 
               */
          }
    }
    cpy[0]=gph[1][0]; cpy[1]=gph[1][1];
    gph[1][0]=gph[r][0]; gph[1][1]=gph[r][1];
    gph[r][0]=cpy[0]; gph[r][1]=cpy[1];
    gph[0][0]=gph[no_edges][0]; gph[0][1]=gph[no_edges][1];
} /* simplify_process2 */


void   Print_Ntks_500(ntks_struct *new_ntks, short rets,  FILE *out_file){
       int graph[MAXSIZE][2];
       short j;
       ntks_struct *ptr;
       short leaf_map[NO_LEAVES];
       short  ntk_p[NO_NODES][THREE];
       short no_nodes, no_edges;
       int no_ntks=0;
  
       

        /* printf("=== Insice print_ntks_500\n"); */
        ptr=new_ntks; 

       while (ptr!=NULL) {
            no_ntks +=1;
            no_edges=ptr->no_edges;
            for (j=0; j<no_edges; j++) {
               graph[j+1][0]=ptr->lt_ends[j]; graph[j+1][1]=ptr->rt_ends[j];
            }
            
            Simplify_Ntks_Aug22(graph, no_edges);

            for (j=1; j<=no_edges; j++) {
              fprintf(out_file, "%d %d\n", graph[j][0], graph[j][1]);
            }
          ptr=ptr->next;
       }

      printf("\n\n    %d ntks  with %d rets are found\n", no_ntks, rets);

} /* print_nks_1 */




ntks_struct *Convert_Zip_to_Gset(g_set *new_ntks){
    g_set *ptr;
    ntks_struct *ptr1, *out_ntks;
    short i;
    int u, v;
    short edge_no;

   /*  printf("inside convert\n"); */
    ptr=new_ntks;
    out_ntks=NULL;
    while (ptr!=NULL) {
       ptr1=(ntks_struct *)malloc(sizeof(ntks_struct));
       edge_no=ptr->no_edges;
       ptr1->lt_ends=(int *)malloc((edge_no+1)*sizeof(int));
       ptr1->rt_ends=(int *)malloc((edge_no+1)*sizeof(int));
       for (i=0; i<edge_no; i++) {
          (ptr1->lt_ends)[i]=(int)(ptr->left_ends)[i]; 
          (ptr1->rt_ends)[i]=(int)(ptr->right_ends)[i]; 
           /* printf("%d %d\n", (ptr->left_ends)[i], (ptr->right_ends)[i]); */
       }
       ptr1->no_edges=edge_no;
       ptr1->next=out_ntks;
       out_ntks=ptr1;
      ptr=ptr->next;
  }
  return out_ntks;
}



/* used */
/* compute the set of leaves below each node of a tree */
void Compute_Paritions(short no_tree_nodes, short no_lf, short sorted_nodes[], 
 char *node_clusters[],short ch[][THREE]){
   short i, j, k;

   for (i=no_tree_nodes; i>=1; i--) {
          k=ch[sorted_nodes[i]][0];
       if (k==0) { 
           node_clusters[sorted_nodes[i]][sorted_nodes[i]]='1';
       } else {
        /*  sorted_nodes[i] is an internal node */
        for (j=1; j<=no_lf; j++) {
            if (node_clusters[ch[sorted_nodes[i]][1]][j]=='1' ||
                node_clusters[ch[sorted_nodes[i]][2]][j]=='1') 
              node_clusters[sorted_nodes[i]][j]='1'; 
        } 
     }
  } /* for */
}

/* used */
void Update(short no_tree_nodes, short no_lfs, char all_zeros[], short  common_part_no[], 
   char  *T0_clusters[], char *node_clusters[]) {
    int i, j;

    for (i=1; i<=no_tree_nodes; i++) {
       for (j=1; j<=no_tree_nodes; j++) {
          if (strcmp(T0_clusters[i], all_zeros)!=0 &&
              strcmp(T0_clusters[i], node_clusters[j])==0) common_part_no[i] +=1;
       }
    }
} /* Update */

/* used */
/* Cluster1 is contained in cluster2, now remove the leaves in clsuter1 from cluster 2 */ 
void Remove_part(short no_lfs, char cluster2[], char cluster1[]){
    int i;

    for (i=1; i<=no_lfs; i++){ if (cluster1[i]=='1') cluster2[i]='0'; } 
}

/* used */
/* check if cluster2 is contained in cluster 1 */
short  Contain_Cluster(short no_lfs, char *cluster1, char *cluster2){
    int i;

    for (i=1; i<=no_lfs; i++){ if (cluster2[i]=='1' && cluster1[i]=='0') return 0; } 
    return 1;
}


/*
 * construction each part of subtrees 
 *  no_t_nodes: the number of nodes in each input tree
 *  subtree holds the reconstructed subtrees 
 * used
*/

short  Construct_Subtrees(short no_t_nodes, short no_lfs, short subtree[][2], 
   short t_ch[][MAX_DEG], short sorted_nodes[], char *node_clusters[],
   char *node_cluster0[],  short index,  short ComponentGraph[][NO_LEAVES] ){

   short i, j, k, p;
   short ans, ans1;
   short flag;
   

   k=0;
   flag=0;
   for (i=1; i<no_t_nodes; i++) {
    ans=Contain_Cluster(no_lfs, node_cluster0[index], node_clusters[sorted_nodes[i]]);

    ans1=0;
    for (j=1; j<= ComponentGraph[index][0]; j++){
     ans1=ans1+Contain_Cluster(no_lfs, node_clusters[sorted_nodes[i]], 
          node_cluster0[ComponentGraph[index][j]]);
    }

     if (ans==1 && ans1>1) {
       if (flag==0) {
        k=1;
        subtree[1][0]=0; subtree[1][1]=sorted_nodes[i];
        flag=1;
       }
       j=t_ch[sorted_nodes[i]][0];
        /*
        printf("     %d   ans1(%d)  %d\n", sorted_nodes[i], ans1, j);  
         */
       for (p=1; p<=j; p++){
         k=k+1;
         subtree[k][0]=sorted_nodes[i]; subtree[k][1]=t_ch[sorted_nodes[i]][p]; 
         /* printf("k(%d)\n", k); */
       }
     } /* if */ 
   } /* i */
   return k;
}  /* construct subtrees */

void top_sort_Aug24(short *gph[], short no_ntk_nodes, 
 short no_ntk_edges, short sorted_nodes[]){
    short i, j, k;
    short  *in_deg=NULL, *ot_deg=NULL, *roots=NULL;
    /* short  child[NO_NODES][MAX_DEG]; */
   short (*child)[MAX_DEG];
         /* 0 contains the number of children */
    short  r_no, no, u;



    in_deg=(short *)malloc(sizeof(short)*(no_ntk_nodes+1));
    ot_deg=(short *)malloc(sizeof(short)*(no_ntk_nodes+1));
    roots=(short *)malloc(sizeof(short)*(no_ntk_nodes+1));
    child=
       (short (*)[MAX_DEG])malloc(sizeof(short [MAX_DEG])*(no_ntk_nodes+1));
    for (i=0; i<=no_ntk_nodes; i++) { in_deg[i]=0;ot_deg[i]=0; }

    for (i=1; i<=no_ntk_edges; i++) {
        in_deg[gph[1][i]] +=1;
        child[gph[0][i]][ot_deg[gph[0][i]]]=gph[1][i];
            /* kth child is child[i][k] */
        ot_deg[gph[0][i]] +=1;
    }

    r_no=0;
    for  (i=0; i<=no_ntk_nodes; i++) {
       if (in_deg[i]==0 &&  ot_deg[i]==1)
           { roots[r_no]=i; r_no=1; break; }  /* root */
    }

    no=0;
    while (r_no >0) {
      u=roots[r_no-1]; sorted_nodes[no]=u;
      no +=1; r_no=r_no-1;
      for (i=0; i<ot_deg[u]; i++) {
         in_deg[child[u][i]]=in_deg[child[u][i]]-1;
         if (in_deg[child[u][i]]==0) {roots[r_no]=child[u][i]; r_no +=1; }
      }
    } /* while loop */
}

void top_sort_Dec20(short graph[][2], short no_ntk_nodes, 
 short no_ntk_edges, short sorted_nodes[]){
    short i, j, k;
    short  *in_deg=NULL, *ot_deg=NULL, *roots=NULL;
    /* short  child[NO_NODES][MAX_DEG]; */
   short (*child)[MAX_DEG];
         /* 0 contains the number of children */
    short  r_no, no, u;



    in_deg=(short *)malloc(sizeof(short)*(no_ntk_nodes+1));
    ot_deg=(short *)malloc(sizeof(short)*(no_ntk_nodes+1));
    roots=(short *)malloc(sizeof(short)*(no_ntk_nodes+1));
    child=
       (short (*)[MAX_DEG])malloc(sizeof(short [MAX_DEG])*(no_ntk_nodes+1));
    for (i=0; i<=no_ntk_nodes; i++) { in_deg[i]=0;ot_deg[i]=0; }

    for (i=1; i<=no_ntk_edges; i++) {
        in_deg[graph[i][1]] +=1;
        child[graph[i][0]][ot_deg[graph[i][0]]]=graph[i][1];
            /* kth child is child[i][k] */
        ot_deg[graph[i][0]] +=1;
    }

    r_no=0;
    for  (i=0; i<=no_ntk_nodes; i++) {
       if (in_deg[i]==0 &&  ot_deg[i]==1)
           { roots[r_no]=i; r_no=1; break; }  /* root */
    }

    no=0;
    while (r_no >0) {
      u=roots[r_no-1]; sorted_nodes[no]=u;
      no +=1; r_no=r_no-1;
      for (i=0; i<ot_deg[u]; i++) {
         in_deg[child[u][i]]=in_deg[child[u][i]]-1;
         if (in_deg[child[u][i]]==0) {roots[r_no]=child[u][i]; r_no +=1; }
      }
    } /* while loop */
}


/* used for decomposition of input trees */
void Common_Clusters(short no_trees, short Input_trees[][NO_NODES][2], 
    short no_tree_nodes,  short no_tree_edges, short no_lfs,
    char *node_clusters[][NO_NODES], short common_part_no[],
    short sorted_tnodes[][NO_NODES], short tree_ch[][NO_NODES][THREE]){
   
    short tree_p[NO_TREES][NO_NODES]; 
    short  leaf_mapping[NO_LEAVES];
    int i, j;
    char all_zeros[NO_LEAVES];

    for (i=0; i<no_trees; i++) {
     /*  renumber nodes int: root=0; internal nodes are numbered from lf_lo+2; */
     Simplify_Ntks(Input_trees[i], no_tree_edges, no_tree_nodes,  leaf_mapping);
     /* compute child and parent relationship: t_ch and t_p */
     Process_Tree(Input_trees[i], no_tree_edges,  no_tree_nodes,  tree_ch[i], tree_p[i]);
   }


   for (j=0; j<=no_lfs; j++) all_zeros[j]='0'; all_zeros[no_lfs+1]='\0';

   /* process tree by tree */
   for (i=0; i<no_trees; i++) {
      top_sort_Dec20(Input_trees[i], no_tree_nodes, no_tree_edges, sorted_tnodes[i]);


     for (j=0; j<=no_tree_nodes; j++) {
          node_clusters[i][j]=(char *)malloc((no_lfs+2)*sizeof(char));
          strcpy(node_clusters[i][j], all_zeros);
     }
     Compute_Paritions(no_tree_nodes,no_lfs,sorted_tnodes[i],node_clusters[i],tree_ch[i]);

     if (i==0) {
        for (j=0; j<=no_tree_nodes; j++){ common_part_no[j]=1; }
     } else {
        Update(no_tree_nodes, no_lfs, all_zeros, common_part_no, node_clusters[0], 
           node_clusters[i]);
     } 
  } /* i loop */

} /* Common_Partition */



void Print_subtrees(short tree[][2], short subtree_size){
   short i;
   
   for (i=1; i<=subtree_size; i++) {
    printf("%d %d\n", tree[i][0], tree[i][1]);
   }
}

/* used for decomponstion of input trees */
/*  ComponentGph[i][0]=no of nodes below i whose clusters are common */
/*  componentGph[i][j] contains the nodex whose cluster is common and contained in 
 *  the partitiion of node i.
 */
void  Compute_ComponentGph(short no_trees, short no_tree_nodes, short no_lfs, 
   short component_graph[][NO_LEAVES], short common_part_no[], short sorted_list[], 
   char *node_clusters[]){

    short i, j, k;
    char  aux_clusters[NO_LEAVES];


   for (i=1; i<no_tree_nodes; i++) {
     k=0;
     if (common_part_no[sorted_list[i]]==no_trees){
        strcpy(aux_clusters, node_clusters[sorted_list[i]]);
        for (j=i+1; j<no_tree_nodes; j++) {
          if ( common_part_no[sorted_list[j]]==no_trees
           && Contain_Cluster(no_lfs, aux_clusters, node_clusters[sorted_list[j]])==1){
              /* the j-th node is below the i-th  */

             Remove_part(no_lfs, aux_clusters, node_clusters[sorted_list[j]]);
            k=k+1;
            component_graph[sorted_list[i]][k]=sorted_list[j];
          } /* if */
        }/* for */
     }/* if */
     component_graph[sorted_list[i]][0]=k;
   } /* i loop */

} /* Compute_ComponentGph */


void  Replacement_Leaves(int *g_left, int *g_right,  short leaf_mapping[], 
   short no_ntk_edges, short no_lfs, short comp_ind, int *copy_left, 
   int *copy_right){
     short i, k;
     int  h0, h1, h;

     /*  printf("inside Replace_Lfs\n"); */

      for (i=0; i<no_ntk_edges; i++) {
         if (g_left[i] > no_lfs)  
           copy_left[i]=g_left[i]+1000*comp_ind;
         else copy_left[i]=g_left[i];
         if (g_right[i] > no_lfs)  copy_right[i]=g_right[i]+1000*comp_ind;
         else  copy_right[i]=g_right[i]; 
      }

     for (i=0; i<no_ntk_edges; i++) {
         h0=copy_left[i]; h1=copy_right[i];
         if (h0==0) {  h=h1; }
         if (1<=h1 && h1<=no_lfs) {
           copy_left[i]=h0;
           copy_right[i]=leaf_mapping[h1];
         }
     }

     for (i=0; i<no_ntk_edges; i++) {
         h0=copy_left[i]; h1=copy_right[i];
         if (h0==h) {
             copy_left[i]=leaf_mapping[0];
             copy_right[i]=h1;
         }
         if (h1==h) { 
            copy_left[i]=h0;
            copy_right[i]=leaf_mapping[0];
         }

        /*  printf("%d %d\n", copy_left[i], copy_right[i]); */
     }
} /* Replacement_Leaves */


ntks_struct *Ntks_Replacement(ntks_struct *ntks, short leaf_mapping[], 
  short no_lfs, short comp_ind){
   ntks_struct *ptr,  *new_list, *new;
   /* int copy_edge[MAXSIZE]; */
   short i;
   short no_ntk_edges;

   /*
   */

   ptr=ntks; new_list=NULL; new=NULL;
   while (ptr!= NULL) {
      no_ntk_edges=ptr->no_edges;

     if (new_list==NULL) { 
          new_list=(ntks_struct *)malloc(sizeof(ntks_struct));
          new_list->no_edges=no_ntk_edges;
          new_list->lt_ends=(int *)malloc(sizeof(int)*no_ntk_edges);
          new_list->rt_ends=(int *)malloc(sizeof(int)*no_ntk_edges);
     Replacement_Leaves(ptr->lt_ends, ptr->rt_ends,  leaf_mapping,  
     no_ntk_edges, no_lfs, comp_ind, new_list->lt_ends, new_list->rt_ends);
   
          new_list->next=NULL; 
          new=new_list;
      } else { 
          new->next=(ntks_struct *)malloc(sizeof(ntks_struct));; 
          (new->next)->no_edges=no_ntk_edges;
          (new->next)->lt_ends=(int *)malloc(sizeof(int)*no_ntk_edges);
          (new->next)->rt_ends=(int *)malloc(sizeof(int)*no_ntk_edges);
     Replacement_Leaves(ptr->lt_ends, ptr->rt_ends,  leaf_mapping,  
     no_ntk_edges, no_lfs, comp_ind, (new->next)->lt_ends, (new->next)->rt_ends);
          new->next->next=NULL;
          new=new->next;
     }
     ptr=ptr->next;
   }

   return new_list;
} /* Ntks_Replacement */


/* attached edges in new_ntks_cp in each ntk of the final_ntks_list */
ntks_struct *Extend_Final_Ntks(ntks_struct *final_ntks_list,  
     ntks_struct *new_ntks_cp, short node_ind){
   short i, no_edges, current_size;
   ntks_struct *ptr, *new, *copy_ptr, *partial_ptr;
   

   /* printf("Inside Extend_Final_Ntks\n"); */
   ptr=NULL;

   if (final_ntks_list==NULL) {
      partial_ptr=new_ntks_cp;
      while (partial_ptr !=NULL) {
         new=(ntks_struct *)malloc(sizeof(ntks_struct));
         no_edges=partial_ptr->no_edges;  
          new->lt_ends=(int *)malloc(sizeof(int)*no_edges);
          new->rt_ends=(int *)malloc(sizeof(int)*no_edges);
        if (node_ind >1) {
          for (i=0; i<no_edges-1; i++){  
            /* new->encoded_edge[i]=(partial_ptr->encoded_edge[i+1]); */
             new->lt_ends[i]=partial_ptr->lt_ends[i+1]; 
             new->rt_ends[i]=partial_ptr->rt_ends[i+1]; 
          }
          new->no_edges=no_edges - 1;
        } else if (node_ind ==1)  {
          for (i=0; i<no_edges; i++){ 
              /* new->encoded_edge[i]=(partial_ptr->encoded_edge[i]); */
             new->lt_ends[i]=partial_ptr->lt_ends[i]; 
             new->rt_ends[i]=partial_ptr->rt_ends[i]; 
          }
           new->no_edges=no_edges;
        }
       /* new->no_edges=no_edges; */
       new->next=ptr; ptr=new; partial_ptr=partial_ptr->next;
     }
  } else {
      copy_ptr=final_ntks_list; 
       
     while (copy_ptr!=NULL){
       current_size=copy_ptr->no_edges;
       partial_ptr=new_ntks_cp;
       while (partial_ptr !=NULL) {
         new=(ntks_struct *)malloc(sizeof(ntks_struct));
         no_edges=partial_ptr->no_edges;
         new->lt_ends=(int *)malloc(sizeof(int)*(current_size+no_edges));
         new->rt_ends=(int *)malloc(sizeof(int)*(current_size+no_edges));
         for (i=0; i<current_size; i++){ 
             new->lt_ends[i]=copy_ptr->lt_ends[i]; 
             new->rt_ends[i]=copy_ptr->rt_ends[i]; 
         }
        if (node_ind ==1){
          for (i=0; i<no_edges; i++) {
             new->lt_ends[current_size+i]=partial_ptr->lt_ends[i]; 
             new->rt_ends[current_size+i]=partial_ptr->rt_ends[i]; 
          }
          new->no_edges=current_size+no_edges;
        } else if (node_ind>1) {
          for (i=0; i<no_edges-1; i++){ 
            new->lt_ends[current_size+i]=partial_ptr->lt_ends[i+1];
            new->rt_ends[current_size+i]=partial_ptr->rt_ends[i+1];
          }
          new->no_edges=current_size+no_edges-1;
        }

          new->next=ptr; ptr=new; partial_ptr=partial_ptr->next;
      }
      copy_ptr=copy_ptr->next;
    } /* while */
 } /* non_empty final_list */

   return ptr;
}


void Adjust_map(short tree[][2], short ind,  short no_lfs, 
  short leaf_mapping[][NO_LEAVES+1], char *node_clusters[][NO_NODES]){
   
     short map[NO_LEAVES]; /*  map[j]=i meaning jth leaf in tree[j]=the jth leaf of tree[0] */
     short map_copy[NO_LEAVES];
     short i, j, k; 
     short n0, ni;
     short  no_tree_edges;

    no_tree_edges=2*no_lfs-1;

     for (i=1; i<=no_lfs; i++) {
         ni=leaf_mapping[ind][i];
         for  (j=1; j<=no_lfs; j++) {
            n0=leaf_mapping[0][j]; 
            if (strcmp(node_clusters[ind][ni], node_clusters[0][n0])==0) map[i]=j;
         }
     }

    for (i=1; i<=no_tree_edges; i++) { k=tree[i][1]; if (k<=no_lfs) tree[i][1]=map[k];   }
    for (i=1; i<=no_lfs; i++) { map_copy[i]=leaf_mapping[ind][i]; }
    for (i=1; i<=no_lfs; i++) { leaf_mapping[ind][map[i]]=map_copy[i]; }
}

   
/* used */
void   Decompose_Trees_Into_Componets(short no_trees, short I_trees[][NO_NODES][2],
   short I_node_no, short I_edge_no, short I_lf_no, char *node_clusters[][NO_NODES], 
   short sorted_tnodes[][NO_NODES], short I_tree_ch[][NO_NODES][THREE], 
    short ComponentGraph[][NO_LEAVES]){

   short common_part_no[NO_NODES];

   Common_Clusters(no_trees, I_trees, I_node_no, I_edge_no,
      I_lf_no, node_clusters, common_part_no, sorted_tnodes, I_tree_ch);

   Compute_ComponentGph(no_trees, I_node_no, I_lf_no, ComponentGraph, common_part_no, 
       sorted_tnodes[0], node_clusters[0]);

}

/* node_index > 1 mean it is the topest component */
ntks_struct *Combine_Solutions_On_Parts(ntks_struct *final_ntks_list, 
  g_set *comp_ntks,  short no_lfs,  
    short leaf_mapping[], short comp_ind,  short node_ind){

    ntks_struct *comp_ntks_cp, *comp_ntks_cpcp;

    comp_ntks_cp=Convert_Zip_to_Gset(comp_ntks);
    comp_ntks_cpcp=Ntks_Replacement(comp_ntks_cp, leaf_mapping, no_lfs,
      comp_ind); 
    return Extend_Final_Ntks(final_ntks_list,  
      comp_ntks_cpcp, node_ind);
  }

void Generate_Permut(char *zh, short pos,  short new_pos, char *current_zh){
      strcpy(current_zh, zh); 
      current_zh[new_pos]=zh[pos]; 
      current_zh[pos]=zh[new_pos];
}

short Min(short a, short b){
   if (a <b) return a; else return b;
}


void Compute_Incidence_July30(short permute[], short inci[NO_NODES],  
   short  pa[], short lfs_no){
    short i, j, k;
    short taxa, parent;

    
    j=2*lfs_no+1;
    for (i=1; i<NO_NODES; i++)  inci[i]=-1;
    inci[0]=permute[1];
    
    for (i=1; i<=lfs_no; i++) {
      taxa=permute[i];
      while (pa[taxa]!=-1) {
        parent=pa[taxa];
        if (inci[parent]==-1) { inci[parent]=permute[i]; taxa=parent;
        } else { inci[parent]=permute[i]; break; }
      } /* while */
   } /* taxa */
}


int ConverteIndex_Aug21(int denom[], short new_pos[], short dim){
    int result;
    short i;

    result=new_pos[0];
    for (i=1; i<dim; i++) { result=result+new_pos[i]*denom[i-1]; }
    return result;
} 


void Int_to_Table_Pos_Aug21(int num,  short dim,  int denom[], short pos[]){
     int r, q;
     short i;

     q=num;
     for (i=1; i<dim; i++) {
       pos[dim-i]=(short)(q/denom[dim-1-i]); 
       q=q%denom[dim-1-i]; 
     } 
     pos[0]=(short)q; 
}


short Inverse(short a, short m,  short pi[]){
   short i;

   for (i=1; i<=m; i++) { if (pi[i]==a) return i; }
}


void Compute_Pos(short new_pos[], short pos[], short appear_index[], short dim){
    short i, num;

    for (i=0; i<dim; i++) new_pos[i]=pos[i];
    num=appear_index[0];
    for (i=1; i<=num; i++) new_pos[appear_index[i]]=pos[appear_index[i]]-1;
}


short Compute_Align_Aug29(short l_pos, short r_pos, short lfs_no, short lens[], 
   short *reduced_seqs[], short **res){

    char *score;
   short *taxa_appear[NO_LEAVES]={NULL};
   short  j, k, taxa, dim;
   short new_pos[NO_LEAVES], pos[NO_LEAVES];
   short no;
   short s, s1;
   int  i, L, num, cell, opt_cell, opt_taxa;
   short seq_len;
   struct inform {
     char taxa;   /* taxa at the current position */
     int   prvs_cell; /* back track purpose */
   } *cond;
   int *denom=NULL;
   /* short mm; */


  dim=r_pos-l_pos+1;
  denom=(int *)malloc(sizeof(int)*dim);
  for (j=0; j<dim; j++) {
     if (j==0) denom[j]=(lens[j+l_pos]+1);
     else denom[j]=denom[j-1]*(lens[j+l_pos]+1);
   }

  score=(char *)malloc((denom[dim-1]+1)*sizeof(char));
  cond=(struct inform *)malloc((denom[dim-1]+1)*sizeof(struct inform));
  for (j=0; j<=lfs_no; j++)
     taxa_appear[j]=(short *)malloc(sizeof(short)*(dim+1));

  score[0]='0';
  L=denom[dim-1]-1;
for (num=1; num<=L; num++) {
      Int_to_Table_Pos_Aug21(num, dim, denom,  pos);

     /* compute the different taxa in that position */
     /* no is the number of different taxa */
       for (j=1; j<=lfs_no; j++) { taxa_appear[j][0]=0; }
       for (j=0; j<dim; j++){
        if (pos[j]>0){
          taxa=reduced_seqs[j+l_pos][pos[j]];
          k=taxa_appear[taxa][0];
          taxa_appear[taxa][0]= k+1;
          taxa_appear[taxa][k+1]=j;
        }
       } /* which taxa in which tree's decomposition */

     /* taxa=1; */
     s=dim*lfs_no+1;
     for (i=1; i<=lfs_no; i++) {
       if (taxa_appear[i][0] >0) {
         Compute_Pos(new_pos, pos, taxa_appear[i], dim);
         cell=ConverteIndex_Aug21(denom, new_pos, dim);
          s1=(short)(score[cell]-'0');
         if ( s>s1  ) { s=s1;  opt_cell=cell; opt_taxa=i; }
       }
     }/* i loop */
     score[num]='0'+(s+1);
     cond[num].taxa='0'+opt_taxa; cond[num].prvs_cell=opt_cell;

  }/* num loop */
  
  num=L; seq_len=0;
  *res=(short *)malloc(sizeof(short)*(2+score[L]-'0'));
  while (num !=0) {
      (*res)[seq_len]=(short)(cond[num].taxa-'0');
      num= cond[num].prvs_cell;
      seq_len++;
  }

  free(score); score=NULL;  free(cond); cond==NULL;
  for (j=0; j<=lfs_no; j++)
      { free(taxa_appear[j]); taxa_appear[j]=NULL;}
  return seq_len;
} /* Comput_Align_Augh29 */

/* this version only find the len of supersequecen, not itself */
short Compute_Score_HighDim(short dim,  short tip, short lfs_no, 
  short *incidence[][NO_LEAVES], short non0_trs1[],  short **align_res){
   char *score;
   short *taxa_appear[NO_LEAVES];
   short new_pos[NO_LEAVES];
   short pos[NO_LEAVES];
   short no;
   short s, s1;
   short  i,  j, k, taxa;
   int  num, cell, L, S, d;
   short len, seq_len, l_pos, r_pos;
   /* short goodness; */
   short *res=NULL;
   short *reduced_seqs[NO_TREES]={NULL};
   short lens[NO_TREES];
   int *denom, max_cells;
   /* MAX INT=2,147,483,647 */


  for (i=0; i<dim; i++){ 
       lens[i]=incidence[non0_trs1[i]][tip][0];
       reduced_seqs[i]=(short *)malloc((lens[i]+1)*sizeof(short));
       for (j=0; j<=lens[i]; j++){ 
           reduced_seqs[i][j]=incidence[non0_trs1[i]][tip][j];
       }
  }

   /* max_cells=100,000,000; */
  max_cells=MEMORY_LIMIT;
   S=1; l_pos=0;
   for (i=0; i<=dim; i++) {
     if (i==dim){
       r_pos=i-1;
       if (r_pos==l_pos) {
         /* len=reduced_seqs[r_pos][0]; */
         len=lens[i];
         if (*align_res!=NULL) free(*align_res);
         *align_res=(short *)malloc(sizeof(short)*(len));
          for (j=1; j<=len; j++) (*align_res)[j-1]=reduced_seqs[r_pos][len-j+1];
          /* format reduced seqs: 0 for size; other positions for  letters  */
       } else {
         len=Compute_Align_Aug29(l_pos, r_pos, lfs_no, lens, reduced_seqs,
            align_res);
        /*
         printf("--Compute_Align_Aug29, l(%d) r(%d) len(%d)\n", l_pos, r_pos, len);
         for (j=0; j<len; j++) { printf("%d ", (*align_res)[j]); }
         printf("\n", dim);
         */
       }
     return len;
     } else if  (max_cells/S < lens[i]+1){ 
       /* Align them  */
       r_pos=i-1;
       len=Compute_Align_Aug29(l_pos, r_pos, lfs_no, lens, reduced_seqs, 
            &res);
       free(reduced_seqs[r_pos]);
       lens[r_pos]=len;
       reduced_seqs[r_pos]=(short *)malloc((len+1)*sizeof(short));
       reduced_seqs[r_pos][0]=len;
       for (j=1; j<=len; j++) { 
            reduced_seqs[r_pos][j]=(res)[len-j];
       }
        S=len*(lens[i]+1); l_pos=r_pos;
     } else  { 
        S=S*(lens[i]+1);
     }
   } /* for */

       for (j=0; j<dim; j++){ 
         if (reduced_seqs[j]!=NULL) free(reduced_seqs[j]);
       }


} /* computer score */

 short  Compute_Score(short dim, short lfs_no,  short lens[],
  short *reduced_seqs[],  short **seq){
   char *score;
   short *taxa_appear[NO_LEAVES];
   short  j, k, taxa;
   /*
   short taxa_appear[NO_LEAVES][NO_LEAVES];
    */
   short new_pos[NO_LEAVES], pos[NO_LEAVES];
   short no;
   short s, s1;
   int  i, L, num, cell, opt_cell, opt_taxa;
   short seq_len;
   struct inform {
     char taxa;   /* taxa at the current position */
     int   prvs_cell; /* back track purpose */
   } *cond;
   int *denom=NULL;
   /* short mm; */
   
  denom=(int *)malloc(sizeof(int)*(dim+1));
  for (j=0; j<dim; j++) {
     if (j==0) denom[j]=(lens[j]+1);
     else denom[j]=denom[j-1]*(lens[j]+1);
   }

  score=(char *)malloc((denom[dim-1]+1)*sizeof(char));
  cond=(struct inform *)malloc((denom[dim-1]+1)*sizeof(struct inform));
  for (j=0; j<NO_LEAVES; j++)
     taxa_appear[j]=(short *)malloc(sizeof(short)*(dim+1));

  score[0]='0';
  L=denom[dim-1]-1;
  for (num=1; num<=L; num++) {
      Int_to_Table_Pos_Aug21(num, dim, denom,  pos);

     /* compute the different taxa in that position */
     /* no is the number of different taxa */
       for (j=1; j<=lfs_no; j++) { taxa_appear[j][0]=0; }
       for (j=0; j<dim; j++){
        if (pos[j]>0){  
          taxa=reduced_seqs[j][pos[j]];
          k=taxa_appear[taxa][0];
          taxa_appear[taxa][0]= k+1;
          taxa_appear[taxa][k+1]=j;
        }
      } /* which taxa in which tree's decomposition */

     /* taxa=1; */
     s=dim*lfs_no+1;
     for (i=1; i<=lfs_no; i++) {
       if (taxa_appear[i][0] >0) {
         Compute_Pos(new_pos, pos, taxa_appear[i], dim);
         cell=ConverteIndex_Aug21(denom, new_pos, dim);
          s1=(short)(score[cell]-'0');
         if ( s>s1  ) { s=s1;  opt_cell=cell; opt_taxa=i; }
       }
     }/* i loop */

     score[num]='0'+(s+1); 
     cond[num].taxa='0'+opt_taxa; cond[num].prvs_cell=opt_cell;

  }/* num loop */


  num=L; seq_len=0;
  if (*seq!=NULL) { free(*seq); *seq=NULL; }
  *seq=(short *)malloc(sizeof(short)*(1+score[L]-'0'));
  /* printf("----I I I here seq_len: %d \n", 1+score[L]-'0'); */
  while (num !=0) {
   /* printf("---%d here %d \n", num, seq_len); */
      (*seq)[seq_len]=(short)(cond[num].taxa-'0'); 
      num= cond[num].prvs_cell; 
   seq_len++;
  }

  free(score); score=NULL;  free(cond); cond==NULL;
  for (j=0; j<NO_LEAVES; j++) 
      { free(taxa_appear[j]); taxa_appear[j]=NULL;}
  return seq_len;
} /* computer score */


short Sum(short pos[], short dim){
  short i, sum;

   sum=0;
  for (i=0; i<dim; i++) { sum=sum+pos[i]; } 
  return sum;
} 

/* two sequence alignment */
short Align_TwoIntSeq(short seq1[], short seq2[]){
   short i, j, k;
   short score[NO_NODES][NO_NODES];
   short len1, len2, len;

   len1=seq1[0]; len2=seq2[0];
   for (i=0; i<=len1; i++) {
     for (j=0; j<=len2; j++) {
        if (i==0) score[i][j]=j;
        else if (j==0) score[i][j]=i;
        else if (seq1[i]==seq2[j]) score[i][j]=1+score[i-1][j-1];
        else score[i][j]=1+Min(score[i][j-1], score[i-1][j]);
     }
   } /* i loop */

   /*
   i=len1; j=len2; len=0;

   while (i>0 && j>0) {
     if (seq1[i]==seq2[j]) { seq[len]=seq1[i]; j--; i--; }
     else if (score[i][j-1]< score[i-1][j]) {
       seq[len]=seq2[j]; j--;
     } else { seq[len]=seq1[i]; i--; }
     len++;
   }

   while (i>0) { seq[len]=seq1[i]; i--; len++;}
   while (j>0) { seq[len]=seq2[j]; j--; len++;}
   */
   return score[len1][len2];
}


short Equal_Containment(short list1[], short list2[]){
   short k;
   
   if (list1[0]<list2[0]) k=list2[0];
   else k=list1[0];

   if (Align_TwoIntSeq(list1, list2)==k) return 1;
   else return 0;

}

/* we assume the fist list is shorter */
short Is_Subsequence(short list1[], short list2[]){
  short i, s1, s2;
  short k;

  s1=list1[0]; s2=list2[0];
  k=1;
  for (i=1; i<=s2; i++) {
    if (list2[i]==list1[k]) k++;
    if (k==s1+1) return 1;
  }
  
  return 0;

}

short Equal_Containment_Aug21(short list1[], short list2[]){
   short i, j, k;

  if (list1[0]<list2[0]) { 
        return Is_Subsequence(list1, list2);
  } else {
        return Is_Subsequence(list2, list1);
  }

}

/* return o if it is not contained in any sequence after that */
short Contain(short *incidence[][NO_LEAVES], short non0_trs[],
   short  id,  short  dim, short taxa){

  short  i;

  for (i=id; i<dim; i++) {
    if (Equal_Containment_Aug21(incidence[non0_trs[id-1]][taxa],
          incidence[non0_trs[i]][taxa])==1) 
     return 1;
  } 
 return 0;

} 


void Convert_Into_Array(short *inci_str[], 
      short taxa_inci[], short lf_no, short p[]) {

    short i, k; 
    short pa, j;
    short  str[NO_LEAVES];

    /* printf("       inside Convert\n"); */
    inci_str[0]=NULL;
    for (i=1; i<=lf_no; i++) {
      pa=p[i]; j=0;
      while (taxa_inci[pa]!=i) {
        str[j]=taxa_inci[pa]; j++; pa=p[pa];
      }
      if (inci_str[i]!=NULL) {free(inci_str[i]); inci_str[i]=NULL;} 
      inci_str[i]=(short *)malloc(sizeof(short)*(j+1));
      inci_str[i][0]=j;
      for (k=1; k<=j; k++) inci_str[i][k]=str[k-1];
    }
    /* printf("          out convert\n"); */
}


short Align2_July8(short no_trees, short taxa, 
  short *incidence[][NO_LEAVES], short lfs_no, short **align_res){
  short i, j, m;
  short len[NO_TREES], non0_trs[NO_TREES], pos[NO_TREES];
  short non0_trs1[NO_TREES];
  short min;
  short dim, index;
  short seq_len;
 /*
  short *reduced_seqs[NO_TREES]={NULL};
  */

  dim=0; 
  for (i=0; i<no_trees; i++) {
     len[i]=(-1)*incidence[i][taxa][0];
     if (incidence[i][taxa][0]>0) { non0_trs[dim]=i; dim++; }
  } /* non0_trs contains all trees where the first components >0)*/


  if (dim==0) { *align_res=NULL;  return 0; }
  else if (dim==1) {
    seq_len=incidence[non0_trs[0]][taxa][0];
    *align_res=(short *)malloc(sizeof(short)*(seq_len)); 
    for (j=1; j<=seq_len; j++) 
      (*align_res)[j-1]=incidence[non0_trs[0]][taxa][seq_len-j+1];  
    return seq_len; 
  }

   MergeSort(non0_trs, 0, dim-1, len);

  /*
  for (i=0; i<dim; i++) {
    printf("In align2_july18: i(%d), %d   %d\n", i, non0_trs[i], 
      len[non0_trs[i]]*(-1));
  }
  printf(" \n\n");
  */

  i=0; 
  for (j=0; j<dim; j++) {
    if (Contain(incidence, non0_trs, j+1,  dim, taxa)==0){
       non0_trs1[i]=non0_trs[j];
       i++;
    }
  }
 /* non0_trs1[i]=non0_trs[dim-1]; */
  dim=i;

  if (dim==1) {
      seq_len=len[non0_trs1[0]]*(-1);
      *align_res=(short *)malloc(sizeof(short)*(seq_len+1)); 
      for (j=1; j<=seq_len; j++) 
        (*align_res)[j-1]=incidence[non0_trs1[0]][taxa][seq_len-j+1];  
      return seq_len;
  }

   /*
  for (i=0; i<dim; i++) {
    printf("align:  i(%d) %d %d: ", i, non0_trs1[i], 
      len[non0_trs1[i]]*(-1));
    
    for (j=1; j<=len[non0_trs1[i]]*(-1); j++){ 
       printf("%d ", incidence[non0_trs1[i]][taxa][j]);
    }
    printf("\n");
  }
  */
 

  seq_len=Compute_Score_HighDim(dim,taxa,lfs_no,incidence, non0_trs1,align_res);
  /*
  printf("   Align2_July out\n");
  */

  return seq_len;
  
} /* Align2_July8 */

/* return -1 indicate  no solution: just find shortest supersequence */
short Align2(short no_trees, short taxa, 
  short *incidence[][NO_LEAVES], short lfs_no, short **seq, short p[]){
  short i, j,  m;
  short sizes[NO_TREES], non0_trs[NO_TREES], pos[NO_TREES];
  short len[NO_TREES], non0_trs1[NO_TREES];
  short base, min;
  short dim, index;
  short seq_len;
  short *reduced_seqs[NO_TREES]={NULL};


  dim=0; 
  for (i=0; i<no_trees; i++) {
     len[i]=(-1)*incidence[i][taxa][0];
     if (incidence[i][taxa][0]>0) { non0_trs[dim]=i; dim++; }
  } 

  if (dim==0) return 0;
  else if (dim==1) {
    seq_len=incidence[non0_trs[0]][taxa][0]; 
    if (*seq!=NULL) { free(*seq); *seq=NULL;}
    *seq=(short *)malloc(sizeof(short)*(seq_len));
    for (i=1; i<=seq_len; i++) 
        (*seq)[i-1]=incidence[non0_trs[0]][taxa][seq_len-i+1];
    return seq_len;
  }

  MergeSort(non0_trs, 0, dim-1, len);

  i=0; 
  for (j=0; j<dim; j++) {
    if (Contain(incidence, non0_trs, j+1,  dim, taxa)==0){
       non0_trs1[i]=non0_trs[j];
       i++;
    }
  }
  dim=i;

  for (i=0; i<dim; i++){ 
     sizes[i]=incidence[non0_trs1[i]][taxa][0];
     if (reduced_seqs[i]!=NULL) { 
         free(reduced_seqs[i]); reduced_seqs[i]=NULL;
     }
     reduced_seqs[i]=(short *)malloc((sizes[i]+1)*sizeof(short));
     for (j=0; j<=sizes[i]; j++)
           reduced_seqs[i][j]=incidence[non0_trs1[i]][taxa][j];
  }


  /* printf("\n\n max(%d) dim(%d), L: %d\n", max, dim,  L); */

  seq_len=Compute_Score(dim, lfs_no,  sizes, reduced_seqs, seq);
  
  for (i=0; i<dim; i++){ free(reduced_seqs[i]); }
  return seq_len;
  
} /* Align2 */
  

short Pairwise_Align_July8(char *seq1, char *seq2){
   short i, j;
   short **score;
   short len1, len2, len;

   len1=strlen(seq1); len2=strlen(seq2);
   score=(short **)malloc(sizeof(char *)*(len1+1));
   for (i=0; i<=len1; i++) score[i]=(short *)malloc(sizeof(short)*(len2+1));

   for (i=0; i<=len1; i++) {
     for (j=0; j<=len2; j++) {
        if (i==0) {score[i][j]=j;
        }
        else if (j==0) { score[i][j]=i;
        }
        else if (seq1[i-1]== seq2[j-1]) {
             score[i][j]=1+score[i-1][j-1];
        }
        else { score[i][j]=1+Min(score[i][j-1], score[i-1][j]);
        }
     }
   } /* i loop */

   j=score[len1][len2];
   for (i=0; i<=len1; i++) free(score[i]);
   free(score);
   return j;
 
}



/* two sequence alignment */
short Pairwise_Align(char *seq1, char *seq2, char seq_str[]){
   short i, j, k;
   short **score;
   short len1, len2, len;
   char *seq;
   char t;

   /* printf("inside: PW align\n"); */
   len1=strlen(seq1); len2=strlen(seq2);
   score=(short **)malloc(sizeof(char *)*(len1+1));
   for (i=0; i<=len1; i++) score[i]=(short *)malloc(sizeof(short)*(len2+1));

   for (i=0; i<=len1; i++) {
     for (j=0; j<=len2; j++) {
        if (i==0) score[i][j]=j;
        else if (j==0) score[i][j]=i;
        else if (seq1[i-1]==seq2[j-1]) score[i][j]=1+score[i-1][j-1];
        else score[i][j]=1+Min(score[i][j-1], score[i-1][j]);
     }
   } /* i loop */
 
   /* printf("score %d\n", score[len1][len2]); */

   i=len1; j=len2; len=0;

   seq=(char *)malloc(score[len1][len2]+1);

   while (i>0 && j>0) {
     if (seq1[i-1]==seq2[j-1]) { seq[len]=seq1[i-1]; j--; i--; }
     else if (score[i][j-1]< score[i-1][j]) { 
       seq[len]=seq2[j-1]; j--; 
     } else { seq[len]=seq1[i-1]; i--; }
     len++;
   }

   while (i>0) { seq[len]=seq1[i-1]; i--; len++;}
   while (j>0) { seq[len]=seq2[j-1]; j--; len++;}

   seq[len]='\0'; 
   for (i=0; i<len; i++) { seq_str[i]=seq[len-1-i]; }
   seq_str[len]='\0';
  /*  printf("out: PW align\n"); */
   for (i=0; i<=len1; i++) free(score[i]);
   free(score);
   return len;
}

short Is_Contain_Str(char *str,  char *sub_str) {
    short i, s1, s2;
  short k;

  s2=strlen(str);
  s1=strlen(sub_str);
  k=0;
  for (i=0; i<s2; i++) { if (sub_str[k]==str[i]) k++; if (k==s1) return 1; }

  return 0;

}

short  Equality_Contain_str(char *seqs_str[], short k,  char *str){ 
     short i;

     for (i=0; i<=k; i++) {
        /* if (strcmp(seqs_str[i], str)==0  ) return 1; */ 
        if (Is_Contain_Str(seqs_str[i], str) ==1) return 1;
     }
     return 0;
}

short Align3(short no_trees, short taxa, char *incidence[][NO_LEAVES], 
short lfs_no, short **seq){
  short i, j, k, m;
  /* short sizes[NO_TREES], non0_trs[NO_TREES]; */
  /* short  pos[NO_TREES]; */
  short min;
  short dim,  dim1;
  short seq_len;
  short q_flag, coord[TWO], co[TWO];
  char  *seqs_str[NO_TREES]={NULL};
  char  *seqs_str1[NO_TREES]={NULL};
  char  *align_result=NULL;
  short *len=NULL, *array=NULL;
  short **scores=NULL;


  dim1=0; 
  for (i=0; i<no_trees; i++) {
     k=strlen(incidence[i][taxa]);
     if (k>0) {
       seqs_str1[dim1]=(char *)malloc(k+1); 
       strcpy(seqs_str1[dim1], incidence[i][taxa]);
       dim1++;
     }
  } /* non0_trs contains all trees where the first components >0)*/

  /* in the following case dela dim>-2 */
  if (dim1==0) {  
     if (seq!=NULL) *seq=NULL;
     return 0;
  }
  if (dim1==1) { 
      seq_len=strlen(seqs_str1[0]);
      if (seq!=NULL) { 
         if (*seq!=NULL) {free(*seq); *seq=NULL;}
         *seq=(short *)malloc(sizeof(short)*seq_len);
         for (j=0; j<seq_len; j++) (*seq)[j]=1+(seqs_str1[0][seq_len-1-j]-'0');
      }
      return seq_len;
  }


  len=(short *)malloc(dim1*sizeof(short));
  array=(short *)malloc(dim1*sizeof(short));
  for (i=0; i<dim1; i++) { 
     len[i]=strlen(seqs_str1[i]);
       array[i]=i;
  }
  MergeSort(array, 0, dim1-1, len);

  i=0; j=1;
  seqs_str[0]=(char *)malloc(strlen(seqs_str1[array[0]])+1);
  strcpy(seqs_str[0], seqs_str1[array[0]]);
  
  /* remove identical copies */
  while (j<dim1) {
    if (Equality_Contain_str(seqs_str, i,  seqs_str1[array[j]])==0){ 
      i=i+1;
      seqs_str[i]=(char *)malloc(strlen(seqs_str1[array[j]])+1);
      strcpy(seqs_str[i], seqs_str1[array[j]]);
    }
    j=j+1;
  }
  dim=i+1;
  free(len);
  free(array);
  for (i=0; i<dim1; i++) free(seqs_str1[i]);
 
  if (dim==1) { 
     seq_len=strlen(seqs_str[0]); 
     if (seq!=NULL) { 
         if (*seq!=NULL) {free(*seq); *seq=NULL;}
         *seq=(short *)malloc(sizeof(short)*seq_len);
         for (j=0; j<seq_len; j++) (*seq)[j]=1+(seqs_str1[0][seq_len-1-j]-'0');
     }
     return seq_len;
 }


  scores=(short **)malloc(sizeof(short *)*(dim+1));
  for (i=0; i<dim; i++) 
  scores[i]=(short *)malloc(sizeof(short)*(dim+1));


  co[0]=-1; co[1]=-1;
  while (dim >2) {
   q_flag=0;
   for (i=0; i<dim; i++) { 
     for (j=i+1; j<dim; j++) {
       /* just pairwise align for the min len */
        if (co[0]==-1 || i==co[0] || j==co[0]) 
            scores[i][j]=Pairwise_Align_July8(seqs_str[i], seqs_str[j]); 
         else if (i==co[1]){ 
            scores[i][j]=scores[j][dim];
         } else if (j==co[1]) scores[i][j]=scores[i][dim]; 
       if (q_flag==0) { min=scores[i][j]; coord[0]=i; coord[1]=j; q_flag=1;
       } else { 
         if (min > scores[i][j]) { min=scores[i][j]; coord[0]=i; coord[1]=j; }
       }
     }
   }
  align_result=(char *)malloc(strlen(seqs_str[coord[0]])+strlen(seqs_str[coord[1]])+1);
   seq_len=Pairwise_Align(seqs_str[coord[0]], seqs_str[coord[1]], align_result);  
   free(seqs_str[coord[0]]);
   seqs_str[coord[0]]=(char *)malloc(strlen(align_result)+1);
   strcpy(seqs_str[coord[0]], align_result);
  free(align_result);

   if (coord[1] !=dim-1) {
      free(seqs_str[coord[1]]);
      seqs_str[coord[1]]=(char *)malloc(strlen(seqs_str[dim-1])+1);
      strcpy(seqs_str[coord[1]], seqs_str[dim-1]);
   }
   co[0]=coord[0]; co[1]=coord[1];
   dim--;
  } /* while */

  /* finding the supersequence */
  if (seq!=NULL) {
    if (*seq!=NULL) { free(*seq); *seq=NULL;}
     *seq=(short *)malloc(sizeof(short)*(strlen(seqs_str[0])+strlen(seqs_str[1])));
     align_result=(char *)malloc(strlen(seqs_str[0])+strlen(seqs_str[1])+1);
     seq_len=Pairwise_Align(seqs_str[0], seqs_str[1], align_result);  
     for (j=0; j<seq_len; j++) (*seq)[j]=1+(align_result[seq_len-1-j]-'0'); 
  } else seq_len=Pairwise_Align_July8(seqs_str[0], seqs_str[1]);
  for (i=0; i<dim; i++) { free(seqs_str[i]); free(scores[i]);}
  free(scores);
  return seq_len;
} /* Align3 */
  

/* in the following construction, the leaves are 1 to lfs_no;
 * the potential ret are lfs_no+2,..., 2lfs_no. tree nodes are
 *  2lfs_no+1... 
 */
/* construct tree-child ntks and return the no. of edges */
short Compute_NTK(short no_trees, short lfs_no, short ret_bound,
 short permut[], short  *inci[][NO_LEAVES], 
  char  *incidence[][NO_LEAVES], short *gph1[2], short *no_nodes){
     short a, b, i, j, k, m, u, v;
     short seq_len, seq_len1;
     short no_edges, node;
     short ret;
     short *seq=NULL;
     short *indeg=NULL;
     short *parent=NULL;
     short *gph[TWO];
     short no_deg1_nodes;
     short **align_res=NULL;

   /*
     printf("in: comput ntk, %d \n", 3*lfs_no+2*ret_bound);
    */

     gph[0]=(short *)
        malloc(sizeof(short)*(3*lfs_no+2*ret_bound));
     gph[1]=(short *)
        malloc(sizeof(short)*(3*lfs_no+2*ret_bound));

    no_edges=0; node=2*lfs_no+1;

   if (align_res!=NULL) { free(align_res); align_res=NULL;}
         align_res=(short **)malloc(sizeof(char *));


    for (i=1; i<lfs_no; i++) {
        m=permut[i];
       
       seq_len=Align2_July8(no_trees, m, inci,  lfs_no, align_res);
       if (no_trees > 2){ 
          seq_len1=Align3(no_trees, m, incidence,  lfs_no, &seq);
          if (seq_len1 < seq_len) { 
            seq_len=seq_len1;
            free(*align_res);
            *align_res=(short *)malloc(sizeof(short)*seq_len);
            for (j=0; j<seq_len; j++)  (*align_res)[j]=seq[j];
          }
        }
       /*
       printf("-----%d---m(%d)-----lf(%d)- %d--\n", i, m, lfs_no, seq_len);
       if (seq_len>0) { 
            for (j=0; j<seq_len; j++) { printf("%d ",  (*align_res)[j]);
         }
         printf("\n\n");
       } else {
       }
       */
       

       if (i==1) gph[0][no_edges]=0; else gph[0][no_edges]=lfs_no+i; 
       /* no_edges++; */

       for (j=0; j<seq_len; j++) {
         gph[1][no_edges]=node; 
         gph[0][no_edges+1]=node; 
         gph[1][no_edges+1]=lfs_no+Inverse((*align_res)[j], lfs_no,  permut); 
         no_edges=no_edges+2;
         gph[0][no_edges]=node; node=node+1;
       }
       gph[1][no_edges]=m;

       no_edges++;
    } /* for  */
       gph[0][no_edges]=lfs_no+lfs_no;
       gph[1][no_edges]=permut[lfs_no];
       no_edges++;

    if (align_res!=NULL) { free(align_res); align_res=NULL;}

    
    parent=(short *)malloc(sizeof(short)*(2*lfs_no+2));
    indeg=(short *)malloc(sizeof(short)*(2*lfs_no+2));
    /* remove degree-two nodes */ 
    a=lfs_no+1; b=2*lfs_no+1;
    for (i=a; i<b; i++) { indeg[i]=0; } 
    for (j=0; j<no_edges; j++) {
       m=gph[1][j]; 
       if (m>=a && m< b) {
           k=indeg[m]+1; indeg[m]=k; 
           if (k==1) parent[m]=gph[0][j];
       }
    } 

   no_deg1_nodes=0;
   for (i=a; i<b; i++) { if (indeg[i]==1) no_deg1_nodes++; }
   *no_nodes=node-1-no_deg1_nodes;  

    k=0;
    for (j=0; j<no_edges; j++) {
        u=gph[0][j]; v=gph[1][j];
        if (v <a || v>=b){ 
          if (u>=a &&  u<b && indeg[u]==1) {
            gph1[0][k]=parent[u]; gph1[1][k]=v; k++;
          } else { gph1[0][k]=u; gph1[1][k]=v; k++;}
        } else  {
          if (indeg[v]==1) {
          } else {
             gph1[0][k]=u; gph1[1][k]=v; k++;
          }
        } 
    }
    gph1[0][k]=gph1[0][0]; gph1[1][k]=gph1[1][0];

    /*
    printf("---------------zzz--------- new edges %d\n", k);
    for  (j=0; j<k; j++) {
       printf("%d %d\n", gph1[0][j], gph1[1][j]);
    }
    printf("\n\n");
    */
    
    free(indeg); free(parent);
    return k;
}


/* return ret number */
short Search_minNTKS_Round1(short no_trees, short lfs_no, short permut[], 
char *incidence[][NO_LEAVES], char *zh, short pos){
     short i;
     short  seq_len, no_rets, no_ret_edges;

    no_ret_edges=0; 
    for (i=1; i<=pos; i++) {
       /* t=permut[i]; */
       seq_len=Align3(no_trees, permut[i], incidence,  lfs_no, NULL);
       no_ret_edges=no_ret_edges+seq_len;
    }

    no_rets=no_ret_edges-pos;
    return no_rets;
    
} /* search min_ntks */

/* new version we use the exact align score */
short Search_minNTKS_Round2(short no_trees, short lfs_no, short permut[], 
short  *inci[][NO_LEAVES], char *incidence[][NO_LEAVES], short pos){
     short i, j;
     short  seq_len, seq_len1;
     short no_rets, no_ret_edges;
     short **align_res=NULL;
     short *seq=NULL;

    
    no_ret_edges=0; 
         if (align_res!=NULL) { free(align_res); align_res=NULL;}
         align_res=(short **)malloc(sizeof(char *));
    for (i=1; i<pos; i++) {
       /* t=permut[i]; */
       seq_len=Align2_July8(no_trees, permut[i], inci,  lfs_no, align_res);
       if (no_trees >2 ){
          seq_len1=Align3(no_trees, permut[i], incidence,  lfs_no, &seq);
          if (seq_len1 < seq_len) seq_len=seq_len1;
       }
       
       if (*align_res!=NULL) {free(*align_res); *align_res=NULL;}
       no_ret_edges=no_ret_edges+seq_len;
    }

    no_rets=no_ret_edges-(pos-1);
    free(align_res);
    return no_rets;
    
} /* search min_ntks */


short Is_Isomorphism(short map[], short r_edges1[][2], short t_edges1[][2],
    short r_edges2[][2], short t_edges2[][2], short array_lf1[], short nodes_no,
    short ret_no, short lf_no){

   short i, j, k, c;
   short current_edge[2];
   short  tree_eg_no, ret_edge_no, used;
   short  flag;

   tree_eg_no=ret_no+2*lf_no-1;
   map[0]=-1;
   for (i=1; i<=lf_no; i++) { map[i]=i; }
   for (i=lf_no+1; i<=nodes_no; i++) { map[i]=-1; }

   c=lf_no;  /* in the start, only leaves are mapped */
   used=0;
   flag=0;
  while (c<nodes_no) {
     /* flag = 1+flag; if (flag > 50) exit(100); */
     for (j=0; j<tree_eg_no - used; j++){
         if (map[t_edges1[j][0]]==-1 && map[t_edges1[j][1]]!=-1){
             k=0;
             while (t_edges2[k][1]!=map[t_edges1[j][1]] && k<tree_eg_no) k++;
             if (k==tree_eg_no) {
                 return 0; }
             map[t_edges1[j][0]]=t_edges2[k][0];
             c=c+1;
             if (j<tree_eg_no-used-1) {
                current_edge[0]=t_edges1[j][0]; current_edge[1]=t_edges1[j][1];
                t_edges1[j][0]=t_edges1[tree_eg_no-used-1][0];
                t_edges1[j][1]=t_edges1[tree_eg_no-used-1][1];
                t_edges1[tree_eg_no-used-1][0]=current_edge[0];
                t_edges1[tree_eg_no-used-1][1]=current_edge[1];
                used= used+1;
             }
             break;
         }
     }/* for */
   } /* c-while loop */

      /* check time */
   ret_edge_no=2*ret_no;
   for (i=0; i<ret_edge_no; i++) {
      if (check_edge(map[r_edges1[i][0]], map[r_edges1[i][1]], r_edges2, ret_edge_no)==0) {
        return 0;
      }
   }
   for (i=0; i<tree_eg_no; i++) {
      if (check_edge(map[t_edges1[i][0]], map[t_edges1[i][1]], t_edges2,
          tree_eg_no)==0) {
         return 0;
       }
   }
   return 1;


} /* Is_Isomorphism */

short Is_Tree_Node(short node, short p[][NO_NODES]){
     if (p[node][0]==1) return 1; else return 0;
}

/* we assume parents are sorted from smallest to the larger */
short Is_Equal(short array1[], short array2[]){
   short i, j, k1, k2;

   k1=array1[0]; k2=array2[0];
   if (k1!=k2) return 0;
   else {
     for (i=1; i<=k1; i++) {
        if (array1[i]!=array2[i]) return 0;
     }
   }
   return 1;
}

short Is_Isomorphism_Dec1(short sorted1[], short p1[][NO_NODES],
  short p2[][NO_NODES], short no_nodes, short no_rets, short no_lfs){

   short i; 
   short map[NO_NODES];
   short xxx, node, image, node2;
   short parent1, parent2;

   map[0]=0;
   for (i=1; i<=no_lfs; i++) { map[i]=i;}
   for (i=no_lfs+1; i<=no_nodes; i++) { map[i]=0;}

    xxx=2*no_lfs;
   for (i=no_nodes-1; i>0; i--) {
       node=sorted1[i];
      if ( Is_Tree_Node(node, p1)==1) { /* tree node */
         image=map[node];
         /* if ( image  > xxx &&  image %2 ==0)  return 0; */
         parent1=p1[node][1];
         parent2=p2[image][1];

         if (map[parent1]==0){
           if ( Is_Tree_Node(parent1, p1)==1){ /* parent is trnode*/
                if  ( Is_Tree_Node(parent2,p2)==1) {
                    map[parent1]=parent2;
                 }  else return 0;
           } else {
               if  (Is_Tree_Node(parent2,p2)==1) return 0;
               else  map[parent1]=parent2;
           }
        } else {
          if (map[parent1]!=parent2) return 0;
        }
     } /* if */
   } /* for */

   for (i=no_nodes-1; i>0; i--) {
       node=sorted1[i];
       node2=map[node];
       if (Is_Tree_Node(node, p1)==0){
          if (Is_Equal(p1[node], p2[node2])==0) return 0;
       }
   }
   return 1;
} /* Dec1 */


short  Check_Iso_v2(short *gph[], g_set *ntks, short edge_no, 
   short node_no, short rets_no,  short lfs_no, short pa[][NO_NODES]){
   short  map[NO_NODES];
   g_set *ntk_ptr;
   short i, j,k, ans;
   short pa2[NO_NODES][NO_NODES];
   short *sorted_nodes=NULL;


      /*
      printf("inside Check_Iso_v2: edge_no(%d) node_node(%d) lfs_no(%d)\n", 
      edge_no, node_no, lfs_no);  
      */
      sorted_nodes=(short *)malloc(sizeof(short)*(node_no+1));
      top_sort_Aug24(gph,  node_no, edge_no, sorted_nodes);
      /*
      for (i=0; i<node_no; i++) printf("%d ", sorted_nodes[i]);
      printf("\n");
       */

     ntk_ptr=ntks;
     while (ntk_ptr!=NULL) {
       if (node_no==ntk_ptr->no_nodes){
        for (i=1; i<node_no+1; i++) {
          k=(ntk_ptr->parent)[i][0];
          for (j=0; j<=k; j++) {
              pa2[i][j]=(ntk_ptr->parent)[i][j];  
          }
        }
        ans=Is_Isomorphism_Dec1(sorted_nodes, pa,  pa2, 
           node_no, rets_no, lfs_no);
        if (ans==1) { free(sorted_nodes); return 1; }
      } else {
      }
      ntk_ptr=ntk_ptr->next;
   }
   free(sorted_nodes);
   return 0;
}

g_set *Insert_v2(g_set *ntk_list, short *gph[], short edge_no, short node_no,
  short ret_no, short lfs_no, char *zh){
   g_set *ptr=NULL;
   short i,k,j, end;
   /*
   short sorted_nodes[NO_NODES];
    */
   short ans;
   short pa[NO_NODES][NO_NODES];


    /*
    printf("Insert graph:  xxx edge_no(%d), node_no(%d) ret_no(%d)\n", 
          edge_no, node_no, ret_no);
     */

       Compute_parents(gph, edge_no, node_no,  pa);


    if (ptr!=NULL) {free(ptr); ptr=NULL;}
    if (ntk_list==NULL) {
       ptr=(g_set *)malloc(sizeof(g_set));
       ptr->no_nodes=node_no;
       ptr->no_edges=edge_no;
       ptr->left_ends=(short *)malloc(edge_no*sizeof(short));
       ptr->right_ends=(short *)malloc(edge_no*sizeof(short));
       for (i=0; i<edge_no; i++) {
         (ptr->left_ends)[i]=gph[0][i+1];
         (ptr->right_ends)[i]=gph[1][i+1];
       }
      /* (ptr->parent)=pa; */
       for (i=1; i<node_no+1; i++) { 
           k=pa[i][0];
           (ptr->parent)[i]=(short *)malloc((k+1)*sizeof(short));
           for (j=0; j<=k; j++) (ptr->parent)[i][j]=pa[i][j];
       }
       ptr->next=NULL;
       /* free(pa); */
       return ptr;
     } else {

      if (Check_Iso_v2(gph, ntk_list, edge_no, node_no, ret_no, lfs_no, pa)==0) {
           /* printf("permutation: %s\n", zh); */
           ptr=(g_set *)malloc(sizeof(g_set));
           ptr->no_nodes=node_no;
           ptr->no_edges=edge_no;
           ptr->left_ends=(short *)malloc(edge_no*sizeof(short));
           ptr->right_ends=(short *)malloc(edge_no*sizeof(short));
           for (i=0; i<edge_no; i++){ 
                 (ptr->left_ends)[i]=gph[0][i+1];
                 (ptr->right_ends)[i]=gph[1][i+1];
           }
           for (i=1; i<node_no+1; i++) { 
               k=pa[i][0];
              (ptr->parent)[i]=(short *)malloc((k+1)*sizeof(short));
              for (j=0; j<=k; j++) (ptr->parent)[i][j]=pa[i][j];
            }
            ptr->next=ntk_list;
            return ptr;
      } else {  
       return ntk_list; 
      }
   }/* if */
   
}


void Convert_Str_Format(char *taxa_incidence_str[], short incidence[][NO_LEAVES], short lfs_no){
    short i, j, k;

     taxa_incidence_str[0]=NULL;
      for (i=1; i<=lfs_no; i++) {
      k=incidence[i][0];
      taxa_incidence_str[i]=(char *)malloc((k+1)*sizeof(char));
      for (j=0; j<k; j++) taxa_incidence_str[i][j]='0'+(incidence[i][j+1]-1);
      taxa_incidence_str[i][k]='\0';
   }
}

short Comput_Cutoff(int ntk_counts[], short min, short max, int total){
    short i, j;
    int size;
    int half;

    size=0;
    half=total/2;
    j=0; 
    for (i=min; i<=max; i++) {
       size=size+ntk_counts[i];
       if (size >= half) return i;
   } 
} /* cutt_off */


void Remove_BadPermut(struct node **Q_front, struct node **Q_end, 
  short cut_off,  short min_ret_no, int ctt){

    struct node *Qptr, *f_ptr, *e_ptr, *bad;

     srand( time(0));    
    if (*Q_front==*Q_end) return;  
    Qptr=(*Q_front)->next; f_ptr=*Q_front; e_ptr=*Q_end;
    if (cut_off > min_ret_no){ 
      while (Qptr!=e_ptr) {
        if (Qptr->score <= cut_off) {
          if ( Qptr->score> min_ret_no) {
            if (rand()%ctt ==0) { f_ptr->next=Qptr; f_ptr=Qptr; } 
          } else {
            f_ptr->next=Qptr; f_ptr=Qptr;
          }
        }
        if (f_ptr!=Qptr) {
          bad=Qptr;
        }
        Qptr=Qptr->next;
       /* if (bad!=NULL) free(bad); */
     }
   } else { 
     while (Qptr!=e_ptr) {
        if (rand()%ctt ==0) { f_ptr->next=Qptr; f_ptr=Qptr; }
        Qptr=Qptr->next;
     }
   }
   if (e_ptr->score > min_ret_no) *Q_end=f_ptr;
}

void Convert_Str_Format_July30(char *inci_str[], short taxa_inci[], short lf_no,
 short p[]) {
    short i; 
    short pa, j;
    char str[NO_LEAVES];

 /*   printf("inside conver leaf_no(%d)\n", lf_no); */

    for (i=1; i<=lf_no; i++) {
      pa=p[i]; j=0;
      while (taxa_inci[pa]!=i) { 
        str[j]='0'+(taxa_inci[pa]-1);
        j++; pa=p[pa];
      }
      
      str[j]='\0';
      inci_str[i]=(char *)malloc(j+1);
      strcpy(inci_str[i], str);
    }
}

/* -------------- */

short  Find_Id(char *taxa_names[], short index,  char *p1){
   short i, j;

   for (i=0; i<index; i++) {
     if (strcmp(taxa_names[i], p1)==0)  return i;
   }
   return -1;
}



void  Replace_Taxa_with_ID_July20(char* str, char* taxa_names[],  short index) {
    char* p = NULL;
    short i,id,k;
    short len = 0;
    char new_str[2000]={0};
    char id_str[100]={0};
    char p1[100]={0}, p2[100]={0};
    short low, high;

    len=strlen(str);
    high=0; id=-1; low=0;
    k=0;
   for (high = 0; high <len; high++) {
        /* printf("       current--%c--\n", str[high]); */
        if (str[high] == ',' || str[high] == ')') {
            if (high != low) {
                memcpy(p1, &str[low], high - low);
                p1[high - low]='\0';
                /* printf("-- yyy---%s--\n", p1); */
                id=Find_Id(taxa_names, index, p1);
                /* printf("            id(%d)\n", id); */
                if (id>=0) {
                  sprintf(id_str, "%d", id + 1);
                } else {
                  sprintf(id_str, "%s", p1);
                }
                  strcat(new_str, id_str);
                  k=strlen(new_str); new_str[k]=str[high];
                  k++;

            } else {
                new_str[k]=str[high]; k++;
            }
            low = high+1;
        } else if (str[high] == '(') {
                new_str[k]=str[high]; k++;
                low = high+1;
        } else {
           /*  new_str[k]=str[high]; k++;  */
        }
   }
   strcpy(str, new_str);
}


/* find how much leaves in the input form */
short Find_Lfs_No(char* c) {
    short i, n;

    n = 0;
    for (i = 0; i < strlen(c); i++) { if (c[i] == ',')  n++; }
    return n + 1;
}




/* number the internal nodes */
void Assign_InternalIDs(char* c) {
    short i, len = 0;
    short n = 0, amount;
    char a[2000];

    amount = Find_Lfs_No(c) + 1;
    /*
 *     for (i = 0; i < strlen(c); i++) {  a[i] = c[i]; }
 *         */
    strcpy(a, c);

    for (i = 0; i < strlen(a); i++) { //insert internal node
        if (a[i] == ')') {
                c[len++] = a[i];
                sprintf(c + len++, "%d", amount++);
                len = strlen(c);
           /* } */
        } else { c[len++] = a[i]; }
    }
    c[len]='\0';
}


/*Data search on strings*/
/* Data_search(nw_str, i + 1, 1, strlen(nw_str), a, tree, c1); */
short  Data_search(char* nw_str, short i, short flag, short l, 
   char* str, short tree[NO_NODES][2], short current_edge) {

    short j = 0, k;
    char taxa_id2[DIGIT_NO] = { 0 }, taxa_id3[DIGIT_NO] = { 0 };
    char a[DIGIT_NO] = { 0 };

    if (flag) {
        while (nw_str[i] >= 48 && nw_str[i] <= 57 && i < l)
        {
          sprintf(taxa_id2, "%d", nw_str[i] - 48); strcat(taxa_id3, taxa_id2);
          i++;
        }
        tree[current_edge][0] = (short) atoi(taxa_id3);
        /* (c1 -> flag3)++; */
    } else {
        while (str[i] >= 48 && str[i] <= 57) { a[j++] = str[i--]; }
        for (k = strlen(a) - 1; k >= 0; k--) {
            sprintf(taxa_id2, "%d", a[k] - 48); strcat(taxa_id3, taxa_id2);
        }
        tree[current_edge - 1][1] = (short) atoi(taxa_id3);
        return i + 1;
    }
}


/*   edge_no=Compute_tree(newick_str, tree, taxa_names, index); */
/*Connect parent and child nodes to generate tree*/
short  Compute_tree(char* nw_str, short tree[NO_NODES][2],  char *taxa_names[], short index) {
    char a[DIGIT_NO] = { 0 };
    short i, j, current_edge;
    char id_str[100];

   /* printf("Comput_tree: %d   %s\n", index, nw_str); */

    j=0;
    current_edge=0;
    for (i = 0; i < strlen(nw_str); i++) {
        a[j++] = nw_str[i];
        if (nw_str[i] == ')')
        {
            while (a[--j] != '(')
            {
       Data_search(nw_str, i + 1, 1, strlen(nw_str), a, tree, current_edge);
                current_edge++;
        j = Data_search(nw_str, j - 1, 0, strlen(nw_str), a, tree, current_edge);

            }
        }
    }
    return  current_edge;
}



short Extract_Taxa(char *newick_str, char *names[]){
    short i;
    char taxa_names[NO_LEAVES][DIGIT_NO];
    short  index, low, high, len;

    len=strlen(newick_str);
    if (newick_str[len - 1] != ')') { newick_str[len - 1] = '\0'; }
    high=0; index=0; low=0;
    for (high = 0; high <len;) {
        if (newick_str[high] == ',' || newick_str[high] == '(' ||
           newick_str[high] == ')') {
            if (high != low) {
                memcpy(taxa_names[index], &newick_str[low], high - low);
                names[index]=(char *)malloc(strlen(taxa_names[index])+1);
                strcpy(names[index], taxa_names[index]);
                index++;
            }
            low = high+1;
        }
        high++;
   }
   return index;
} /* Extract_Taxa */


short  Convert_to_Edges(short flag, char *newick_str, char *taxa_names[],
    short tree[][2], short *lfs_no) {

    short i = 0, index, p=0, len, edge_no;
    char no_internalID_format[2000];
    char id_str[100];



   if (flag==0) {
      index=Extract_Taxa(newick_str, taxa_names);
      /* strcpy(no_internalID_format, newick_str); */
      *lfs_no=index;
   }

   Replace_Taxa_with_ID_July20(newick_str, taxa_names, *lfs_no);
   Assign_InternalIDs(newick_str);
   edge_no=Compute_tree(newick_str, tree, taxa_names, *lfs_no);

   return edge_no;
} /* Convert_to_Edges */



void Add_Rootedge(int edge_no, short  tree[][2]){
   short node_no;
   short *in_deg, *out_deg;
   short i, root;

   in_deg=(short *)malloc((edge_no+2)*sizeof(short));
   out_deg=(short *)malloc((edge_no+2)*sizeof(short));

   for (i=0; i<= edge_no+1; i++) {
     in_deg[i]=0; out_deg[i]=0;
   }

   for (i=0; i< edge_no; i++) {
     in_deg[tree[i][1]] +=1; out_deg[tree[i][0]] +=1;
   }

   for (i=0; i<= edge_no+1; i++) {
     if (in_deg[i]==0 && out_deg[i]==2) {root=i; break;}
   }

   tree[edge_no][0]=tree[0][0]; tree[edge_no][1]=tree[0][1];
   tree[edge_no+1][0]=tree[1][0]; tree[edge_no+1][1]=tree[1][1];
   tree[1][0]=0; tree[1][1]=root;
} /* Add_Rootedge */



void main(int argc, char *argv[]){

FILE *input_file, *tr_input_file, *out_file; 
char *node_clusters[NO_TREES][NO_NODES];
int node1, node2, size_xxx;

short i, j, k, k1;
short no_rets, no_lfs, no_ntk_edges, no_ntk_nodes; /* ntkwork parameters */
short no_trees; 
short subtrees[NO_TREES][NO_NODES][TWO];

short I_trees[NO_TREES][NO_NODES][TWO]; /* keep ith tree in I_trees[i] */
short I_lf_no, I_edge_no, I_node_no;
short I_tree_ch[NO_TREES][NO_NODES][THREE];

short subtr_edge_no, subtr_lf_no;
short leaf_mapping[NO_TREES][NO_LEAVES+1];
short sorted_tnodes[NO_TREES][NO_NODES];
short ComponentGraph[NO_NODES][NO_LEAVES]; 
short node_ind, current_comp;
short current_size;
short t_ch[NO_TREES][NO_NODES][THREE], *t_p[NO_TREES];
short tr_sorted[NO_TREES][NO_NODES];
short no_components;

struct node *Q_front, *Q_end, *Qptr;
struct node  *ptr;
struct node **Second_Q[TWO];
short q_flag, min_ret_no, new_min_ret_no, current_ret_no; 
short min_node_no, ret_val;
short permut[NO_LEAVES+1];
char  *zh, *current_zh;
short *gph[2], no_gph_edges, no_gph_nodes;
short taxa_inci[NO_TREES][NO_NODES];
      /* indicate in the tree decomposition, whihc taxa is right below which
         and the depth of their least common ancestor */
short *align_result[NO_LEAVES], align_flag[NO_LEAVES]; 
           /* first component align result */
char *taxa_incidence_str[NO_TREES][NO_LEAVES];
char *taxa_incidence_str1[NO_TREES][NO_LEAVES];
short *taxa_inci_int[NO_TREES][NO_LEAVES];
int ctt;
short max, kkk, cut_off;
int *ntks_counts=NULL;
short  it_times, pos;
short total_rets=0;
short i00, xxx;
long int  count_ntks;

char input_nw_str[2000];
char *taxa_names[NO_LEAVES];

ntks_struct  *final_ntks_list;
g_set *ntk_list;
/* INT_MAX=2,147,483,647 */



   if (argc>=4) {
      I_lf_no=atoi(argv[2]);
      tr_input_file=fopen(argv[1], "r");
      if (tr_input_file==NULL ) { printf("File openning fail!\n"); exit(10); }
   } else { printf("./a.out <tree_input_file> <no_of_taxa> <ouput_file>\n"); exit(100); }


   /* --------------- start reading input trees --------------------  */
   k=0; /* keep the number of edges  read */ 
   no_trees=0;

   /*
   I_edge_no=2*I_lf_no-1; 
   I_node_no=2*I_lf_no;  
    */
    /* including 0, deg-1 root */

   while (fscanf(tr_input_file, "%s\n", input_nw_str)!=EOF) {
      if (no_trees==0) {
      I_edge_no=Convert_to_Edges(0, input_nw_str, taxa_names, I_trees[no_trees], &I_node_no);
      Add_Rootedge(I_edge_no,  I_trees[no_trees]);
      } else {
       Convert_to_Edges(1, input_nw_str, taxa_names, I_trees[no_trees], &I_node_no);
      Add_Rootedge(I_edge_no,  I_trees[no_trees]);
      }
      no_trees++;
   }
   fclose(tr_input_file);
   /* --------------- finish reading input trees --------------------  */
    
      I_edge_no++;
   I_node_no=I_edge_no+1;   /* including 0, deg-1 root */


  /*
  for (i=0; i<no_trees; i++) {
     for (j=1; j<=I_edge_no; j++) {
       printf("%d %d\n", I_trees[i][j][0], I_trees[i][j][1]);
     }
     printf("\n\n");

  }
   */

   out_file=fopen(argv[3], "w"); fclose(out_file);


  /* decomponse subtrees into many disjoint subset of irreducible subtrees */
   Decompose_Trees_Into_Componets(no_trees, I_trees, I_node_no, I_edge_no,
      I_lf_no, node_clusters, sorted_tnodes, I_tree_ch, ComponentGraph);

   no_components=0;
   for (node_ind=I_node_no-1; node_ind>=1; node_ind--)
       if (ComponentGraph[sorted_tnodes[0][node_ind]][0]>0)  no_components +=1;

   current_comp=0; final_ntks_list=NULL; 

  /* examine  the tree nodes in the 0-th trees one by one */
  for (node_ind=I_node_no-1; node_ind>=1; node_ind--){

   /* identify a common component, then work on it */
   if (ComponentGraph[sorted_tnodes[0][node_ind]][0]>0){
      current_comp +=1;
      printf(" COMPONENT %d\n", current_comp);

   
     for (j=0; j<no_trees; j++){
       /* construct the subtrees rooted at the node of the paritition 
         node_clusters[0][sorted_tnodes[0][node_ind]] */
       subtr_edge_no=Construct_Subtrees(I_node_no, I_lf_no, subtrees[j], I_tree_ch[j], 
           sorted_tnodes[j], node_clusters[j], node_clusters[0],
              sorted_tnodes[0][node_ind], ComponentGraph); 

       
    
       subtr_lf_no=(subtr_edge_no+1)/2;
      
       t_p[j]=(short *)malloc(sizeof(short)*(subtr_edge_no+2));

       Simplify_Ntks_Initialization(subtrees[j], subtr_edge_no,  subtr_edge_no+1,
       leaf_mapping[j], current_comp, no_components);
       /* no_components is used */
 
       if (j>0) Adjust_map(subtrees[j], j, subtr_lf_no, leaf_mapping, node_clusters); 
       /* Compute child array and parent array */
       Process_Tree(subtrees[j], subtr_edge_no,  subtr_edge_no+1,   t_ch[j], t_p[j]);
       /* sort nodes in trees */
       top_sort(t_ch[j], 2*subtr_lf_no, subtr_edge_no, tr_sorted[j], subtr_lf_no);


       
       /*
       printf("after  normalization \n\n");
       Print_subtrees(subtrees[j], subtr_edge_no);
       printf("with the following leaf mapping \n");
       for (k=0; k<=subtr_lf_no; k++) printf("   %d %d\n", k, leaf_mapping[j][k]);
       printf("\n\n");
       */
      

     } /* for j loop */



      zh=(char *)malloc(subtr_lf_no+1); 
      for (j=0; j<subtr_lf_no; j++) zh[j]='0'+j; zh[subtr_lf_no]='\0';
      current_zh= (char*)malloc(subtr_lf_no+1);

      Q0_front=NULL; Q0_end=NULL;  Q1_front=NULL; Q1_end=NULL;
      kkk=0; it_times=subtr_lf_no/2;
      /* Qinsert(&Q0_front, &Q0_end, zh); */

      ntks_counts=(int *)malloc(sizeof(int)*(no_trees)*subtr_lf_no);
     for (kkk=0; kkk < it_times; kkk++ ) {
         if (kkk==0) { 
             Qinsert(&Q0_front, &Q0_end, zh); 
             cut_off=(no_trees-1)*subtr_lf_no;}
         else { Q0_front=Q_front; Q0_end=Q_end; }

         while (Qdelete(&Q0_front, &Q0_end, &ret_val,  zh)==1) {
             if (ret_val <= cut_off) {
                for (k1=2*kkk; k1<subtr_lf_no; k1++) {
                   Generate_Permut(zh, 2*kkk, k1, current_zh);
                   /* printf("-- current_zh: top %s\n", current_zh); */  
                   Qinsert(&Q1_front, &Q1_end, current_zh);
                }
             }
         }


         if (kkk !=it_times-1 || subtr_lf_no%2!=0) {
         while (Qdelete(&Q1_front, &Q1_end, &ret_val,  zh)==1) {
             for (k1=(2*kkk+1); k1<subtr_lf_no; k1++) {
                Generate_Permut(zh, 2*kkk+1, k1, current_zh);
                 /*  printf("-- current_zh: %s\n", current_zh); */  
                Qinsert(&Q0_front, &Q0_end, current_zh);
             }
         }
        } /* if condition */

      
        /*
          printf("----2nd------ kkk(%d)---------\n", kkk); 
         */

        Q_front=NULL; Q_end=NULL;
       if (kkk !=it_times-1 || subtr_lf_no%2!=0) { 
         Q_front=Q0_front; Q_end=Q0_end;
         pos=2*kkk+1;
       } else { pos=2*kkk; 
            Q_front=Q1_front; Q_end=Q1_end;
           /*  Qremove(&Q0_front, &Q0_end); */
       }

       if (Q_front!=Q_end) { 
         printf("      count in the middle pos(%d) xxxx\n", pos); 
         /*
         Count_NNN(Q_front, Q_end); 
          */
       }

       q_flag=0; ctt=0; 
       Qptr=Q_front;
     
       for (j=0; j<(no_trees)*subtr_lf_no; j++) ntks_counts[j]=0;   
       

       while (Qptr != NULL) {
          strcpy(zh, Qptr->item);
          ctt++;
          for (j=1; j<=subtr_lf_no; j++){  permut[j]=1+(zh[j-1]-'0'); }
          for (j=0; j<no_trees; j++) {
             Compute_Incidence_July30(permut, taxa_inci[j],  t_p[j],
               subtr_lf_no);
            Convert_Str_Format_July30(taxa_incidence_str[j], taxa_inci[j], 
                subtr_lf_no, t_p[j]); 
          } /* j loop for subtrees */

          current_ret_no=Search_minNTKS_Round1(no_trees, subtr_lf_no, permut,
          taxa_incidence_str, zh, pos+1);




         Qptr->score=current_ret_no;
         ntks_counts[current_ret_no]=1+ntks_counts[current_ret_no];
         if (q_flag==0) {
           max=current_ret_no; min_ret_no= current_ret_no; q_flag=1;
         } else {
           if (max< current_ret_no) max=current_ret_no;
           if (min_ret_no > current_ret_no) min_ret_no= current_ret_no;
         }
         Qptr=Qptr->next;
       } /* queue loop */
   

       cut_off=Comput_Cutoff(ntks_counts, min_ret_no, max, ctt);
        cut_off=Min(min_ret_no+8, cut_off); 
       if (ctt>CUT_OFF_SIZE) {
          ctt=ctt/(2*CUT_OFF_SIZE); 
          if (ctt==0) ctt=1;
          Remove_BadPermut(&Q_front, &Q_end, cut_off, min_ret_no, ctt);
      }
      /* printf("after remove\n");
      Count_NNN(Q_front, Q_end); 
      */
   } /* kkk loop */
       free(ntks_counts); ntks_counts=NULL;
       

   
    printf("     After sampling search--- min_ret_no(%d)\n", min_ret_no);
     /*
     if (Q_front!=Q_end) {  Count_NNN(Q_front, Q_end); }
      */

   count_ntks=0;
   Qptr=Q_front; new_min_ret_no=min_ret_no;

   q_flag=0;
   while (Qptr!=Q_end) {
      if  (Qptr->score ==min_ret_no) {
        strcpy(zh, Qptr->item);


        for (j=1; j<=subtr_lf_no; j++){  permut[j]=1+(zh[j-1]-'0'); }
        for (j=0; j<no_trees; j++) {
          Compute_Incidence_July30(permut, taxa_inci[j],  t_p[j], subtr_lf_no);
       Convert_Into_Array(taxa_inci_int[j], taxa_inci[j], subtr_lf_no, t_p[j]);
            Convert_Str_Format_July30(taxa_incidence_str[j], taxa_inci[j], 
                subtr_lf_no, t_p[j]); 
        } 


       current_ret_no=Search_minNTKS_Round2(no_trees, subtr_lf_no, permut, 
        taxa_inci_int, taxa_incidence_str, subtr_lf_no);

        if (q_flag==0) { 
          new_min_ret_no =current_ret_no; q_flag=1;
   Second_Q[0]=(struct node  **)malloc(sizeof(struct node *)*(min_ret_no+1));
   Second_Q[1]=(struct node  **)malloc(sizeof(struct node *)*(min_ret_no+1));
   for (j=0; j<=new_min_ret_no; j++) { Second_Q[0][j]=NULL; Second_Q[1][j]=NULL;}
         count_ntks=0;
       }


        if (new_min_ret_no > current_ret_no) {
            new_min_ret_no= current_ret_no; count_ntks=1;
            Qinsert(&(Second_Q[0][new_min_ret_no]), &(Second_Q[1][new_min_ret_no]), zh);
        } else if (new_min_ret_no == current_ret_no) {
          
          count_ntks++;
          Qinsert(&(Second_Q[0][new_min_ret_no]), &(Second_Q[1][new_min_ret_no]), zh);
        }
    } /* if */
     Qptr=Qptr->next;
   } /* Q while */ 
   free(zh);

      

      min_ret_no=new_min_ret_no;
      printf("       the current networks with ret_no(%d) after round 2\n", new_min_ret_no);
      /*
      Count_NNN(Second_Q[0][min_ret_no], Second_Q[1][min_ret_no]); 
       */


     ntk_list=NULL;
     zh=(char *)malloc(subtr_lf_no+1);
     gph[0]=(short *)malloc(sizeof(short )*(3*subtr_lf_no+2*min_ret_no));
     gph[1]=(short *)malloc(sizeof(short )*(3*subtr_lf_no+2*min_ret_no));
      xxx=0;
     while (Qdelete(&(Second_Q[0][min_ret_no]),&(Second_Q[1][min_ret_no]), &ret_val, zh)==1){
        /*  printf("start string  %s\n", zh); */  
       for (j=1; j<=subtr_lf_no; j++)  { permut[j]=1+(zh[j-1]-'0'); }
      for (j=0; j<no_trees; j++) { 
           Compute_Incidence_July30(permut, taxa_inci[j],  t_p[j],
               subtr_lf_no);
           Convert_Into_Array(taxa_inci_int[j], taxa_inci[j], subtr_lf_no, t_p[j]);
            Convert_Str_Format_July30(taxa_incidence_str[j], taxa_inci[j], 
                subtr_lf_no, t_p[j]); 
      }

      no_gph_edges=Compute_NTK(no_trees, subtr_lf_no, min_ret_no, permut, 
             taxa_inci_int, taxa_incidence_str,  gph,  &no_ntk_nodes);
      Simplify_Ntks_July22(gph, no_gph_edges);
      ntk_list=Insert_v2(ntk_list, gph, no_gph_edges, no_ntk_nodes,
           no_gph_edges-no_ntk_nodes+1, subtr_lf_no, zh);
    } /* while */
  
    free(gph[0]); gph[0]=NULL;
    free(gph[1]); gph[1]=NULL;
    total_rets=total_rets+no_gph_edges-no_ntk_nodes+1;


      

    final_ntks_list=Combine_Solutions_On_Parts(final_ntks_list, 
             ntk_list, subtr_lf_no, leaf_mapping[0], current_comp, node_ind);

     free(Second_Q[0]); Second_Q[0]=NULL; free(Second_Q[1]); Second_Q[1]=NULL;
   } /* component */
 }    /* node */

 out_file=fopen(argv[3], "a"); 
 if (out_file==NULL) {  printf("cannot open out_file\n"); exit(100); }

 /* print out multiple networks, whose size information can be foudn in 
 *  Progress_Log */
 Print_Ntks_500(final_ntks_list, total_rets,  out_file);
 for (j=0; j<no_trees; j++) free(t_p[j]);
 close(out_file);

} /* end main */
