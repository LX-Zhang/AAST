/*    Copyright@Louxin Zhang  National University of Singapore
   This program takes a set of (weighted, unweighted) rooted binary phylogenetic
   trees and compute tree-child networks (acyclic rooted networks with a set
   of leaves satisfying some constraints)
   that contain every input tree and have the minimum  number of reticulations,
   i.e., no networks with less reticulations contain all the trees.
   Suhc a network has alos the smallest number of nodes among those containing
   all input trees.

   The input trees are listed row by row. Each row contains a tree in Newick format, such as
      ((1,(3,2)7)8)

   is a tree on taxe 1, 2, 3. Here, we assume the root 0 has out-degree 1.

  Compiling command:  gcc -O3 ALTS.c  -o NTK_Finder
  Run command: ./NTK_Finder <tree_input_file> <no_leaves> <output_file>

   Remarks:  1. The networks output in <ouputfile> have the same number of
            nodes. But the network may have different number of edges becasue they can
       have non-binary redticulation nodes, which have 2 or more incoming edges.

      2. The output networks are listed one by one, separated by a line containing
      "\\". each network is listed by edges, one in each row.

     3. IN current setting, the input tree set can contian at most 100
        binary trees and in each tree, there are at most 128 leaves.
       The program can abort because of memeory usage (1G)

     4. In the revised version, the input trees don't necessarily the same taxa.

*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<time.h>

#define TREE_OR_LEAF 0
#define RETNODE 1
#define T_EDGE 'a'
#define RED_EDGE 'b'
#define T_LEAF -1
#define T_ROOT 2
#define T_OTHER 1

#define MAX_DEG 100
#define MAX_OUT  3
#define NO_TREES 100
#define NO_NODES  300  
 /* this needs to be greater than the number of nodes in the output ntk */ 
#define NO_LEAVES 128
#define MAXSIZE  1024   /* deal with ntks with at most 20 leaves */
#define TWO 2
#define  THREE 3
#define TEST 4   /* test purpose */
#define TAXA_NAME_STR  100 
#define TREE_INPUT_STR 2000

#define CUT_OFF_SIZE 4000
#define MEMORY_LIMIT 500000;



/* encode an edge (x, y) by a int integer: z*100+y to save memory space */
typedef struct graph {
    int *left_ends;    
    int *right_ends;    
    int  *parent[NO_NODES];    /* two parents p1*100+p3 */
    int no_nodes;
    int no_edges;
    struct graph *next;
} g_set;

/* a graph is represented by a linked list of edges (x, y), 
 * which are encoded by an interger using * ENCODE(1000)*x+y */
typedef struct f_gph {
   /* int encoded_edge[MAXSIZE]; */
   int *lt_ends;
   int *rt_ends;
   int no_edges;
   struct f_gph *next;
} ntks_struct;

struct node {
   char *item;
   int  score;
   struct node* next;
};

struct node *Q0_front=NULL, *Q1_front=NULL;
struct node *Q0_end=NULL, *Q1_end=NULL;


void Qinsert_23(struct node **front, struct node **end,  char *zh) {
     struct node *ptr;

  /*  printf("Inside Qinsert\n"); */
   
     ptr = (struct node *) malloc (sizeof(struct node));
   if (ptr == NULL) { 
       printf("\nOVERFLOW\n"); return;
   } else {
     ptr->item = (char *)malloc(strlen(zh)+1);
     strcpy(ptr->item, zh); 
     ptr->score=0; 
     ptr->next =NULL;
    if (*front == NULL) { 
         *front = ptr; *end = ptr;
    } else { 
         (*end) -> next = ptr; (*end) = ptr; 
    }
   }
} /* Qinsert_23 */



int  Qdelete(struct node **front, struct node **end, int  *score,  char *str) {  
    struct node *ptr;


    if (*front == NULL) {  
         return 0;  
    }  else {
        ptr=(*front); strcpy(str, ptr->item);  *score=ptr->score;
        if (*front==*end) { 
             *front=NULL; *end=NULL;
        } else { (*front) = ptr -> next;  }
        free(ptr->item);  free(ptr);
        return 1;
    } 
    return 0;
}   
 

/* Merges two subarrays of arr[].  First subarray is arr[left..m]
  Second subarray is arr[m+1..right]
 */
 void Merge(int A[], int left, int  m, int right, int len[]) {
	int i, j, k;
        int  n1 = m - left + 1;
        int  n2 = right - m;
        int  L[n1], R[n2];

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


/* left is for left index and r is right index of 
 * the sub-array of arr to be sorted */
 void MergeSort(int  A[], int left, int right, int len[]) {
     int m;

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
  while (ptr!=B) { 
      i++; 
      ptr=ptr->next; 
  }
  i++;
  printf("         %d ntks\n", i);
}



/* convert a tree or a ntk that is kept as an array of edges into a
 * child array and parent array 
 * Ch[i] keeps the children of the ith node. Ch[i][0] equals 
 * the number of childrens of
 *  the ith node. each tree node has two children, recorded 
 *  in Ch[i][1] and Ch[i][2]; *  the root and each ret have 1, kept in Ch[i][1].
 */
void Process_Tree_Dec2023(int tree[][TWO], int ch[][MAX_OUT], int no_total_lfs){
    int  t, k;
    int no_nodes, edges_no;

   edges_no=tree[0][0];
   no_nodes=2*no_total_lfs;

   for (t=0; t<=no_nodes;  t++) { ch[t][0]=0; }
   for (t=1; t<=edges_no;  t++) {
      ch[tree[t][0]][0] +=1; k=ch[tree[t][0]][0];
      ch[tree[t][0]][k]=tree[t][1];
   }
} /* process100 */


void Process_Tree_Jan2024_v2(int tree[][TWO], int lf_no, int ch[][MAX_OUT], int pa[]){
    int  i, t, k;
    int max_nodes_id;
    int no_edges;
    int  in_deg[NO_NODES],ot_deg[NO_NODES];
    int  internal[NO_NODES];

   
   no_edges=tree[0][0];
   /* max_nodes_id=2*+2; */

    for (i=0; i<NO_NODES; i++) {  ot_deg[i]=0; ch[i][0]=0;  }
    for (i=1; i<=no_edges; i++) { ot_deg[tree[i][0]] +=1;  }
    k=2; 
    for  (i=0; i<NO_NODES; i++) {
       if (ot_deg[i]==2) {  internal[i]=k; k +=1;  /* tree node */ } 
    }
   
   for (i=1; i<=no_edges; i++) {
      if (ot_deg[tree[i][0]]==2) tree[i][0]= lf_no + internal[tree[i][0]];
      if (ot_deg[tree[i][1]]==2) tree[i][1]= lf_no + internal[tree[i][1]];
  }
  
  /*  for (t=0; t<NO_NODES;  t++) { ch[t][0]=0; pa[t]=-1;} */
   for (t=1; t<=no_edges;  t++) {
      ch[tree[t][0]][0] +=1; k=ch[tree[t][0]][0]; 
      ch[tree[t][0]][k]=tree[t][1];
      pa[tree[t][1]]=tree[t][0];
   }
} /* process100 */

void Process_Tree_Jan2024(int tree[][TWO],  int ch[][MAX_OUT], int pa[]){
    int  t, k;
    int max_nodes_id;
    int edges_no;
   
   edges_no=tree[0][0];
   max_nodes_id=edges_no+2;
  
   for (t=0; t<max_nodes_id;  t++) { ch[t][0]=0; pa[t]=-1;}
   for (t=1; t<=edges_no;  t++) {
      ch[tree[t][0]][0] +=1; k=ch[tree[t][0]][0]; 
      ch[tree[t][0]][k]=tree[t][1];
      pa[tree[t][1]]=tree[t][0];
   }
} /* process100 */


/* change indeg to a pointer */
void Compute_parents(int *gph[TWO], int edges_no, int nodes_no, int *pa[NO_NODES]){
    int  i, m, n, k, min, pos;
    int  *in_deg=NULL;

   /* pa[0][0]=-1; */
   in_deg=(int *)malloc(sizeof(int)*(nodes_no+2));
   for (i=0; i<nodes_no+2; i++) in_deg[i]=0;
   for (i=1; i<=edges_no;  i++) { in_deg[gph[1][i]]+=1; }

    pa[0]=(int *)malloc(sizeof(int)); pa[0][0]=-1;

   for (i=1; i<nodes_no+2; i++){ 
         pa[i]=(int *)malloc((in_deg[i]+1)*sizeof(int)); 
         pa[i][0]=0;
   }
   for (i=1; i<=edges_no;  i++) {
      k=pa[gph[1][i]][0];
      pa[gph[1][i]][0]=1+k;  pa[gph[1][i]][k+1]=gph[0][i];
   }

   free(in_deg);
    
} /* Compute parent array */


int Check_clusters(int *no_lfs_top, char *top_clusters[], char *cluster){
   int no, i;

   no=*no_lfs_top+1;
   for (i=1; i<no; i++) {
     if ( strcmp(top_clusters[i], cluster)==0) return i;
   }
   *no_lfs_top=no;
   top_clusters[no]=(char *)malloc(strlen(cluster)+1);
   strcpy(top_clusters[no], cluster);
   return no;
}

void  Simplify_Tree_Initialization_Jan2024(int gph[][2],  int leaf_map[], 
  int *no_lfs_top, char *top_clusters[], char *clusters_array[]){
    int i, j, k, r;
    int  in_deg[NO_NODES],ot_deg[NO_NODES];  
    int  internal[NO_NODES],external[NO_NODES];
    int  v_type[NO_NODES];
    int no_lfs, no_tree_nodes, cpy[2];
    int no_edges;
    int indicator=0;

    no_edges=gph[0][0];
    for (i=0; i<NO_NODES; i++) {  ot_deg[i]=0;  in_deg[i]=0;}
    for (i=1; i<=no_edges; i++) { ot_deg[gph[i][0]] +=1;  in_deg[gph[i][1]]+=1;}
    k=0; j=1;
    for  (i=0; i<NO_NODES; i++) {
       if (  ot_deg[i]==1) { v_type[i]=4;
       }  else if (ot_deg[i]==2) {  
            v_type[i]=3; internal[i]=k; k +=1;  /* tree node */
       } else if (ot_deg[i]==0 &&  in_deg[i]==1) {
               v_type[i]=1;  /* leaf */
               external[i]=j; j +=1;
       } 
    }
    
   no_lfs=j+1; no_tree_nodes=k;
   for (i=1; i<=no_edges; i++) {
          if (v_type[gph[i][0]]==4) { 
             r=i;   gph[i][0]=0; leaf_map[0]=gph[i][1];
           } else if (v_type[gph[i][0]]==3) {} 
              /* gph[i][0]= no_lfs + internal[gph[i][0]]; */

          if (v_type[gph[i][1]]==3){
              /*  gph[i][1]= no_lfs + internal[gph[i][1]]; */
          } else if (v_type[gph[i][1]]==1) { 
   /* if the cluster is foudn, retunr the index; otherwise, attached it and return index */
             indicator=Check_clusters(no_lfs_top, top_clusters, clusters_array[gph[i][1]]); 
             leaf_map[indicator]=gph[i][1];
             gph[i][1]=indicator; 
          }
          /*
          printf("%d %d\n", gph[i][0],  gph[i][1]);
          */
  }

  cpy[0]=gph[1][0]; cpy[1]=gph[1][1];
  gph[1][0]=gph[r][0]; gph[1][1]=gph[r][1];
  gph[r][0]=cpy[0]; gph[r][1]=cpy[1];
  gph[0][0]=no_edges;

} /*end of Simplify_Tree */

/* used */
/* gph is a tree */
void  Simplify_Tree_Initialization(int gph[][2],  int leaf_map[]){
    int i, j, k, r;
    int  in_deg[NO_NODES],ot_deg[NO_NODES];  
    int  internal[NO_NODES],external[NO_NODES];
    int  v_type[NO_NODES];
    int no_lfs, no_trs_rets, cpy[2];
    int no_edges;

    no_edges=gph[0][0];
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
    gph[0][0]=no_edges;

} /* simplify_process2 */


void Simplify_Ntks(int gph[][2], int no_edges, int no_nodes){
    int i, j, k, r;
    int  in_deg[NO_NODES],ot_deg[NO_NODES],internal[NO_NODES];
    int external[NO_NODES], v_type[NO_NODES], cpy[2];
    int no_lfs, no_trs_rets;

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
           } else if (v_type[gph[i][0]]==3 ||v_type[gph[i][0]]==0) 
              /* tree nodes or ret */
              gph[i][0]= no_lfs + internal[gph[i][0]];

          if (v_type[gph[i][1]]==3 || v_type[gph[i][1]]==0)
              gph[i][1]= no_lfs + internal[gph[i][1]];
          else if (v_type[gph[i][1]]==1) { 
               /* gph[i][i] is a leaf, rename it as external[gph[i][1]] */ 
               gph[i][1]=external[gph[i][1]]; 
          }
          /*
          printf("%d %d\n", gph[i][0],  gph[i][1]);
          */
    }

    cpy[0]=gph[1][0]; cpy[1]=gph[1][1]; /* copy the first edge  for swap */
    gph[1][0]=gph[r][0]; gph[1][1]=gph[r][1]; 
          /* move the the edge leaving the root to the first edge */
    gph[r][0]=cpy[0]; gph[r][1]=cpy[1];
              /* move the first edge to the r-th edge */
    gph[0][0]=gph[no_edges][0]; gph[0][1]=gph[no_edges][1]; 
      /* copy the last edge to the 0-th edge */
} /* Simplify_Ntks */



/* In this version,  we consider input trees who may not have the same taxa */
/* rename the  internal nodes from no_total_lfs+2 */
void Simplify_Ntks_Dec2023(int gph[][2], int no_total_lfs){
    int i, j, k, r;
    int  in_deg[NO_NODES],ot_deg[NO_NODES],internal[NO_NODES];
    int external[NO_NODES], v_type[NO_NODES], cpy[2];
    int internal_no1;
    int no_edges, no_nodes;

    no_edges=gph[0][0]; no_nodes=no_edges+1;
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

   internal_no1=no_total_lfs+2;  
   for (i=1; i<=no_edges; i++) {
          if (v_type[gph[i][0]]==4) { 
              r=i;   gph[i][0]=0; 
           } else if (v_type[gph[i][0]]==3 ||v_type[gph[i][0]]==0) 
              /* tree nodes or ret */
              gph[i][0]= internal_no1 + internal[gph[i][0]];

          if (v_type[gph[i][1]]==3 || v_type[gph[i][1]]==0)
              gph[i][1]= internal_no1 + internal[gph[i][1]];
          else if (v_type[gph[i][1]]==1) { 
               /* gph[i][i] is a leaf, rename it as external[gph[i][1]] */ 
               /* in this version, mark off the following sentence, it dones not make
 *                sentence */
               /*  gph[i][1]=external[gph[i][1]];  */
          }
          /*
          printf("%d %d\n", gph[i][0],  gph[i][1]);
          */
    }

    cpy[0]=gph[1][0]; cpy[1]=gph[1][1]; /* copy the first edge  for swap */
    gph[1][0]=gph[r][0]; gph[1][1]=gph[r][1]; 
          /* move the the edge leaving the root to the first edge */
    gph[r][0]=cpy[0]; gph[r][1]=cpy[1];
              /* move the first edge to the r-th edge */
    /*
    gph[0][0]=gph[no_edges][0]; gph[0][1]=gph[no_edges][1]; 
     */
      /* copy the last edge to the 0-th edge */
} /* Simplify_Ntks_Dec23 */

void Print_Ntk(int *gph[], int no_edges){
   int i;


   for (i=1; i<=no_edges; i++) {
     printf("%d %d\n", gph[0][i], gph[1][i]);
   }
    printf("end\n");
}

void  Simplify_Ntks_July22(int *gph[], int no_edges){
    int i, j, k, r;
    int  *in_deg, *ot_deg;
    int *internal;
    int  *v_type, cpy[TWO];
    int no_lfs, no_trs_rets;
    int max_node_id;


    
    max_node_id=0;
    for (i=1; i<=no_edges; i++) { 
            if (gph[1][i]>max_node_id) max_node_id=gph[1][i];
            if (gph[0][i]>max_node_id) max_node_id=gph[0][i];
    }
   

    in_deg=(int*)malloc((max_node_id+2)*sizeof(int));
    ot_deg=(int*)malloc((max_node_id+2)*sizeof(int));
    internal=(int*)malloc((max_node_id+2)*sizeof(int));
    v_type=(int*)malloc((max_node_id+2)*sizeof(int));
    
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
          if (v_type[gph[0][i]]==4) { r=i;   gph[0][i]=0; }
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

    free(in_deg); free(ot_deg); free(internal); free(v_type);
} /* Simplify_Ntks_July22 */


/* compute the number of tree nodes in a network */
int  Simplify_Ntks_Aug22(int gph[][2], int no_edges){
    int i;
    int  j, k, r;
    int  *in_deg, *ot_deg;
    int *internal;
    int  *v_type, cpy[TWO];
    int no_lfs, no_trs_rets;
    int max_node_id;

    max_node_id=0;
    for (i=1; i<=no_edges; i++) { 
            if (gph[i][1]>max_node_id) max_node_id=gph[i][1];
            if (gph[i][0]>max_node_id) max_node_id=gph[i][0];
    }
   
    /* printf("---- max: %d\n", max_node_id); */

    in_deg=(int*)malloc((max_node_id+2)*sizeof(int));
    ot_deg=(int*)malloc((max_node_id+2)*sizeof(int));
    internal=(int*)malloc((max_node_id+2)*sizeof(int));
    v_type=(int*)malloc((max_node_id+2)*sizeof(int));
    
    for (i=0; i<=max_node_id; i++) { in_deg[i]=0; ot_deg[i]=0; }

    for (i=1; i<=no_edges; i++) { 
       in_deg[gph[i][1]] +=1; ot_deg[gph[i][0]] +=1;
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
          if (v_type[gph[i][0]]==4) { r=i;   gph[i][0]=0; 
           }
          if (v_type[gph[i][0]]==3 ||v_type[gph[i][0]]==0) 
              gph[i][0]= no_lfs + internal[gph[i][0]];
          if (v_type[gph[i][1]]==3 || v_type[gph[i][1]]==0)
              gph[i][1]= no_lfs + internal[gph[i][1]];

          if (v_type[gph[i][1]]==1) { }
    }
    cpy[0]=gph[1][0]; cpy[1]=gph[1][1];
    gph[1][0]=gph[r][0]; gph[1][1]=gph[r][1];
    gph[r][0]=cpy[0]; gph[r][1]=cpy[1];
    gph[0][0]=gph[no_edges][0]; gph[0][1]=gph[no_edges][1];

    free(in_deg); free(ot_deg); free(internal); free(v_type);

    return  no_lfs+no_trs_rets-1;
} /* simplify_process2 */


ntks_struct *Convert_Zip_to_Gset(g_set *new_ntks){
    g_set *ptr;
    ntks_struct *ptr1, *out_ntks;
    int i, no;
    /* int u, v; */
    int edge_no;

   /*  printf("inside convert\n"); */
    ptr=new_ntks;
    out_ntks=NULL;
    no=0;
    while (ptr!=NULL) {
       if (no>=2) break;
       ptr1=(ntks_struct *)malloc(sizeof(ntks_struct));
       edge_no=ptr->no_edges;
       ptr1->lt_ends=(int *)malloc((edge_no+1)*sizeof(int));
       ptr1->rt_ends=(int *)malloc((edge_no+1)*sizeof(int));
       for (i=0; i<edge_no; i++) {
          (ptr1->lt_ends)[i]=(ptr->left_ends)[i]; 
          (ptr1->rt_ends)[i]=(ptr->right_ends)[i]; 
       }
       ptr1->no_edges=edge_no;
       ptr1->next=out_ntks;
       out_ntks=ptr1;
      ptr=ptr->next;
      no++;
  }
  return out_ntks;
}



/* used */
/* compute the set of leaves below each node of a tree */
void Compute_Paritions(int no_tree_nodes, int no_all_lf, int sorted_nodes[], 
 char *node_clusters[],int ch[][MAX_OUT]){
   int i, j, k;

   for (i=no_tree_nodes-1; i>=1; i--) {
          k=ch[sorted_nodes[i]][0];
       if (k==0) { 
           /* leaf case */
           node_clusters[sorted_nodes[i]][sorted_nodes[i]]='1';
       } else {
        /*  sorted_nodes[i] is an internal node */
        /* assume degree 2 for internal nodes */
        for (j=1; j<=no_all_lf; j++) {
            if (node_clusters[ch[sorted_nodes[i]][1]][j]=='1' ||
                node_clusters[ch[sorted_nodes[i]][2]][j]=='1') 
              node_clusters[sorted_nodes[i]][j]='1'; 
        } 
     }
  } /* for */
  /*
  printf("in_comput_partition---%d--\n", no_tree_nodes);
  for (i=no_tree_nodes-1; i>=1; i--) printf("--%s--\n", node_clusters[sorted_nodes[i]]);
  printf("\nxxxxx\n");
  */
} /* Compute_Paritions */

/* used */
void Update(int max_node_index,  char *all_zeros, int  common_part_no[], 
   char  *T0_clusters[], char *node_clusters[]) {
    int i, j;

    for (i=1; i<=max_node_index; i++) {
       for (j=1; j<=max_node_index; j++) {
          if (strcmp(T0_clusters[i], all_zeros)!=0 &&
              strcmp(T0_clusters[i], node_clusters[j])==0) common_part_no[i] +=1;
       }
    }
} /* Update */

/* used */
/* Cluster1 is contained in cluster2, now remove the leaves in clsuter1 from cluster 2 */ 
void Remove_part(int no_all_lfs, char cluster2[], char cluster1[]){
    int i;

    for (i=1; i<=no_all_lfs; i++){ if (cluster1[i]=='1') cluster2[i]='0'; } 
}

/* used */
/* check if cluster2 is contained in cluster 1 */
int  Contain_Cluster(int no_all_lfs, char *cluster1, char *cluster2){
    int i;

    for (i=1; i<=no_all_lfs; i++){ 
       if (cluster2[i]=='1' && cluster1[i]=='0') return 0; 
    } 
    return 1;
}


/*
 * construction each part of subtrees 
 *  no_t_nodes: the number of nodes in each input tree
 *  subtree holds the reconstructed subtrees 
 * used
*/

int  Construct_Subtrees(int no_t_nodes, int no_all_lfs, int subtree[][2], 
   int t_ch[][MAX_OUT], int sorted_nodes[], char *node_clusters[],
   char *node_cluster0[],  int index,  int ComponentGraph[][NO_LEAVES] ){

   int i, j, k, p, ans, ans1, flag;
   

   k=0; flag=0;
   for (i=1; i<no_t_nodes; i++) {
    ans=Contain_Cluster(no_all_lfs, node_cluster0[index], node_clusters[sorted_nodes[i]]);

    ans1=0;
    for (j=1; j<= ComponentGraph[index][0]; j++){
     ans1=ans1+Contain_Cluster(no_all_lfs, node_clusters[sorted_nodes[i]], 
          node_cluster0[ComponentGraph[index][j]]);
    }

     if (ans==1 && ans1>1) {
       if (flag==0) {
        k=1;
        subtree[1][0]=0; subtree[1][1]=sorted_nodes[i];
        flag=1;
       }
       j=t_ch[sorted_nodes[i]][0];
       for (p=1; p<=j; p++){
         k=k+1;
         subtree[k][0]=sorted_nodes[i]; subtree[k][1]=t_ch[sorted_nodes[i]][p]; 
         /* printf("k(%d)\n", k); */
       }
     } /* if */ 
   } /* i */
   return k;
}  /* construct subtrees */

int  Construct_Subtrees_Jan2024(int no_all_lfs, int no_nodes,  int subtree[][2], 
   int t_ch[][MAX_OUT], int sorted_nodes[], char *node_clusters[],
   char *node_cluster0[],  int index,  int ComponentGraph[][NO_LEAVES] ){

   int i, j, k, p, ans, ans1, flag;
   

   k=0; flag=0;
   for (i=1; i<no_nodes; i++) {
    ans=Contain_Cluster(no_all_lfs, node_cluster0[index], node_clusters[sorted_nodes[i]]);

    if (ans==1) {
      ans1=0;
      for (j=1; j<= ComponentGraph[index][0]; j++){
        ans1=ans1+Contain_Cluster(no_all_lfs, node_clusters[sorted_nodes[i]], 
          node_cluster0[ComponentGraph[index][j]]);
      }


     if (ans==1 && ans1>1) {
       if (flag==0) {
        k=1;
        subtree[1][0]=0; subtree[1][1]=sorted_nodes[i];
        flag=1;
       }
       j=t_ch[sorted_nodes[i]][0];
       for (p=1; p<=j; p++){
         k=k+1;
         subtree[k][0]=sorted_nodes[i]; subtree[k][1]=t_ch[sorted_nodes[i]][p]; 
         /* printf("k(%d)\n", k); */
       }
     } /* if */ 
    } /* ans */
   } /* i loop */
   subtree[0][0]=k;
   return k;
}  /* construct subtrees */

int Copy_Subtrees_Top_Jan2024(int oritinal_tree[][2], int subtrees[][2]){
  int i, k;
  k=oritinal_tree[0][0];
  for (i=1; i<= k; i++) {
    subtrees[i][0]=oritinal_tree[i][0];
    subtrees[i][1]=oritinal_tree[i][1];
  }
  subtrees[0][0]=k;
  return k;
}

void top_sort_Aug24(int *gph[], int no_ntk_nodes, 
 int no_ntk_edges, int sorted_nodes[]){
    int i, j, k;
    int  *in_deg=NULL, *ot_deg=NULL, *roots=NULL;
   int (*child)[MAX_OUT];
         /* 0 contains the number of children */
    int  r_no, no, u;



    in_deg=(int *)malloc(sizeof(int)*(no_ntk_nodes+1));
    ot_deg=(int *)malloc(sizeof(int)*(no_ntk_nodes+1));
    roots=(int *)malloc(sizeof(int)*(no_ntk_nodes+1));
    child=
       (int (*)[MAX_OUT])malloc(sizeof(int [MAX_OUT])*(no_ntk_nodes+1));
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

     free(in_deg); free(ot_deg); free(roots); free(child);

}

/* revised in Jan 2024 */
void top_sort_Jan2024(int graph[][2], int no_nodes, int no_edges, int sorted_nodes[]){
    int i, j, k;
    int  *in_deg=NULL, *ot_deg=NULL, *indeg_0_nodes=NULL;
   int (*child)[MAX_OUT];
         /* 0 contains the number of children */
    int  r_no, no, u;



    in_deg=(int *)malloc(sizeof(int)*(no_nodes+1));
    ot_deg=(int *)malloc(sizeof(int)*(no_nodes+1));
    indeg_0_nodes=(int *)malloc(sizeof(int)*(no_nodes+1));
    child= (int (*)[MAX_OUT])malloc(sizeof(int [MAX_OUT])*(no_nodes+1));

    for (i=0; i<=no_nodes; i++) { in_deg[i]=0;ot_deg[i]=0; }

    for (i=1; i<=no_edges; i++) {
        in_deg[graph[i][1]] +=1;
        child[graph[i][0]][ot_deg[graph[i][0]]]=graph[i][1];
            /* kth child is child[i][k] */
        ot_deg[graph[i][0]] +=1;
    }

    r_no=0;
    for  (i=0; i<=no_nodes; i++) {
       if (in_deg[i]==0 &&  ot_deg[i]==1)
           { indeg_0_nodes[r_no]=i; r_no=1; break; }  /* root */
    }

    no=0;
    while (r_no >0) {
      u=indeg_0_nodes[r_no-1]; 
      sorted_nodes[no]=u;
      no +=1; r_no=r_no-1;
      for (i=0; i<ot_deg[u]; i++) {
         in_deg[child[u][i]]=in_deg[child[u][i]]-1;
         if (in_deg[child[u][i]]==0) {indeg_0_nodes[r_no]=child[u][i]; r_no +=1; }
      }
    } /* while loop */
     free(in_deg); free(ot_deg); free(indeg_0_nodes); free(child);
}

void top_sort_Dec25(int graph[][2], int no_ntk_nodes, 
 int no_ntk_edges, int no_lfs, int sorted_nodes[]){
    int i, j, k;
    int  *in_deg=NULL, *ot_deg=NULL, *roots=NULL;
   int (*child)[MAX_OUT];
    int  r_no, no, u;

    in_deg=(int *)malloc(sizeof(int)*(no_ntk_nodes+1));
    ot_deg=(int *)malloc(sizeof(int)*(no_ntk_nodes+1));
    roots=(int *)malloc(sizeof(int)*(no_ntk_nodes+1));
    child=
       (int (*)[MAX_OUT])malloc(sizeof(int [MAX_OUT])*(no_ntk_nodes+1));
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
      u=roots[r_no-1]; 
      if (u>no_lfs) { sorted_nodes[no]=u; no +=1;} 
      r_no=r_no-1;
      for (i=0; i<ot_deg[u]; i++) {
         in_deg[child[u][i]]=in_deg[child[u][i]]-1;
         if (in_deg[child[u][i]]==0) {roots[r_no]=child[u][i]; r_no +=1; }
      }
    } /* while loop */
     free(in_deg); free(ot_deg); free(roots); free(child);
}


/* used for decomposition of input trees */
/*  no_tree_nodes= 2*I_lf_no, no_tree_edges= 2*I_lf_no-1 */
/* In the new version, input trees may not have the same taxa. Hence, the no_edges infromation is given the I_tree[][0][0] */
void Common_Clusters(int no_trees, int Input_trees[][NO_NODES][2], 
    int no_all_lfs,
    char *node_clusters[][NO_NODES], int common_part_no[],
    int sorted_tnodes[][NO_NODES], int tree_ch[][NO_NODES][MAX_OUT]){
   
    int  i, j;
    char *all_zeros;
    int max_node_index, no_tree_edges;
    int no_tree_nodes;


    for (i=0; i< no_trees; i++) {
     /*  renumber nodes int: root=0; internal nodes are numbered from lf_lo+2; */
     /* this lookds not necessary as the process in the previous state, Add the root */
     /* Simplify_Ntks_Dec2023(Input_trees[i], no_lfs); */
     /* compute child and parent relationship: t_ch and t_p */
     Process_Tree_Dec2023(Input_trees[i],  tree_ch[i], no_all_lfs);
    /*  i++; */
   }


   all_zeros=(char *)malloc(no_all_lfs+2);
   for (j=0; j<=no_all_lfs; j++) all_zeros[j]='0'; all_zeros[no_all_lfs+1]='\0';

   /* process tree by tree */
    max_node_index=2*no_all_lfs;

   for (i=0; i<no_trees; i++) {
      no_tree_edges=Input_trees[i][0][0]; 
      /* max_node_index=no_all_lfs+(no_tree_edges+1)/2; */
      top_sort_Jan2024(Input_trees[i], max_node_index, no_tree_edges, sorted_tnodes[i]);

      for (j=0; j<=max_node_index; j++) { strcpy(node_clusters[i][j], all_zeros); }

      no_tree_nodes=no_tree_edges+1;
      Compute_Paritions(no_tree_nodes, no_all_lfs,sorted_tnodes[i],node_clusters[i],tree_ch[i]);


     if (i==0) { for (j=0; j<=max_node_index; j++){ common_part_no[j]+=1; }
     } else {
        Update(max_node_index, all_zeros, common_part_no, node_clusters[0], node_clusters[i]);
     } 

  } /* i loop */

} /* Common_Clusters( */




/* used for decomponstion of input trees */
/*  ComponentGph[i][0]=no of nodes below i whose clusters are common */
/*  componentGph[i][j] contains the nodex whose cluster is common and contained in 
 *  the partitiion of node i.
 */
void  Compute_ComponentGph(int no_trees, int no_tree_nodes, int no_all_lfs, 
   int component_graph[][NO_LEAVES], int common_part_no[], int sorted_list[], 
   char *node_clusters[]){

    int i, j, k;
    char  aux_clusters[NO_LEAVES];


   for (i=1; i<no_tree_nodes; i++) {
     k=0;
     if (common_part_no[sorted_list[i]]==no_trees){
        strcpy(aux_clusters, node_clusters[sorted_list[i]]);
        for (j=i+1; j<no_tree_nodes; j++) {
          if ( common_part_no[sorted_list[j]]==no_trees
           && Contain_Cluster(no_all_lfs, aux_clusters, node_clusters[sorted_list[j]])==1){
              /* the j-th node is below the i-th  */

             Remove_part(no_all_lfs, aux_clusters, node_clusters[sorted_list[j]]);
            k=k+1;
            component_graph[sorted_list[i]][k]=sorted_list[j];
          } /* if */
        }/* for */
     }/* if */
     component_graph[sorted_list[i]][0]=k;
   } /* i loop */

   printf("component graph information\n");
   for (i=0; i<no_tree_nodes; i++) {
      printf("--%d--%d--\n", sorted_list[i], component_graph[sorted_list[i]][0]);
      if (component_graph[sorted_list[i]][0]>0) {
       for (j=1; j<= component_graph[sorted_list[i]][0]; j++) 
          printf("xx%d", component_graph[sorted_list[i]][j]);
       printf("\n");
      }
   }

} /* Compute_ComponentGph */



/*  removes all the edges from the orginial trees, and update the number of edges */
void Remove_Subtree_Edges(int orginal_tree[][2],int subtree[][2]){ 
  int i, j, k, u, v;
  int no_edges;

 no_edges=subtree[0][0];
  k=orginal_tree[0][0];
  for (i=2; i<=no_edges; i++) {
     u=subtree[i][0]; v=subtree[i][1];
     for (j=1;  j<=k;  j++) {
       if (u==orginal_tree[j][0] && v==orginal_tree[j][1]){
         if (j==k) { orginal_tree[j][0]=0; orginal_tree[j][1]=0; } else {
            orginal_tree[j][0]=orginal_tree[k][0]; 
            orginal_tree[j][1]=orginal_tree[k][1];
         }
         k=k-1;
         break;
       }
     }
  }
  orginal_tree[0][0]=k;
}



void  Replacement_Leaves(int *g_left, int *g_right,  int leaf_mapping[], 
   int no_ntk_edges, int no_lfs, int comp_ind, int *copy_left, 
   int *copy_right){
     int i, k;
     int  h0, h1, h;


      for (i=0; i<no_ntk_edges; i++) {
         /* printf("%d %d\n", g_left[i], g_right[i]); */
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
             copy_left[i]=leaf_mapping[0]; copy_right[i]=h1;
         }
         if (h1==h) { 
            copy_left[i]=h0; copy_right[i]=leaf_mapping[0];
         }

         /*  printf("%d %d\n", copy_left[i], copy_right[i]);  */
     }
} /* Replacement_Leaves */


ntks_struct *Ntks_Replacement(ntks_struct *ntks, int leaf_mapping[], 
  int no_lfs, int comp_ind){
   ntks_struct *ptr,  *new_list, *new;
   int i;
   int no_ntk_edges;


   ptr=ntks; new_list=NULL; new=NULL;
   while (ptr!= NULL) {
      no_ntk_edges=ptr->no_edges;
      /*  printf("---no_edges(%d)--no_lfs(%d)\n", no_ntk_edges, no_lfs); */

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
          new->next=(ntks_struct *)malloc(sizeof(ntks_struct)); 
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
     ntks_struct *new_ntks_cp, int node_ind){
   int i, no_edges, current_size;
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
             new->lt_ends[i]=partial_ptr->lt_ends[i+1]; 
             new->rt_ends[i]=partial_ptr->rt_ends[i+1]; 
          }
          new->no_edges=no_edges - 1;
          /* this line is modified for unequal case */
        } else if (node_ind ==1 || node_ind ==0)  {
          for (i=0; i<no_edges; i++){ 
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
        /* this line is modified for unequal taxa case */
        if (node_ind ==1 || node_ind ==0){
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


/*
void Adjust_map(int tree[][2], int tree_id, 
  int leaf_mapping[][NO_LEAVES+1], char *node_clusters[][NO_NODES]){
   
     int map[NO_LEAVES];  
     int map_copy[NO_LEAVES];
     int i, j, k; 
     int n0, ni;
     int  no_tree_edges, no_lfs;

     no_tree_edges=tree[0][0];
     no_lfs=(no_tree_edges+1)/2;

     for (i=1; i<=no_lfs; i++) {
         ni=leaf_mapping[tree_id][i]; 
         for  (j=1; j<=no_lfs; j++) {
            n0=leaf_mapping[0][j]; 
            if (strcmp(node_clusters[tree_id][ni], node_clusters[0][n0])==0) map[i]=j;
         }
     }

    for (i=1; i<=no_tree_edges; i++) { k=tree[i][1]; if (k<=no_lfs) tree[i][1]=map[k];   }
    for (i=1; i<=no_lfs; i++) { map_copy[i]=leaf_mapping[tree_id][i]; }
    for (i=1; i<=no_lfs; i++) { leaf_mapping[tree_id][map[i]]=map_copy[i]; }
} 
*/
/* adjust */

   
/* used */
/* remoce the I_node_no, I_edge_no */
void   Decompose_Trees_Into_Componets(int no_trees, int I_trees[][NO_NODES][2],
    int I_lf_no, char *node_clusters[][NO_NODES], 
   int sorted_tnodes[][NO_NODES], int I_tree_ch[][NO_NODES][MAX_OUT], 
    int ComponentGraph[][NO_LEAVES]){

   int common_part_no[NO_NODES];
   int i, no_nodes;

   for (i=0; i<NO_NODES; i++) common_part_no[i]=0;
   Common_Clusters(no_trees, I_trees, 
      I_lf_no, node_clusters, common_part_no, sorted_tnodes, I_tree_ch);
   /* check is done in Jan */

   no_nodes=I_lf_no+(I_trees[0][0][0]+1)/2;
   Compute_ComponentGph(no_trees,  no_nodes,  I_lf_no, ComponentGraph, common_part_no, 
       sorted_tnodes[0], node_clusters[0]);
} /* Decompose_Trees */

/* node_index > 1 mean it is the topest component */
ntks_struct *Combine_Solutions_On_Parts(ntks_struct *final_ntks_list, 
  g_set *comp_ntks,  int no_lfs,  
    int leaf_mapping[], int comp_ind,  int node_ind){

    ntks_struct *comp_ntks_cp, *comp_ntks_cpcp;

    comp_ntks_cp=Convert_Zip_to_Gset(comp_ntks);
 comp_ntks_cpcp=Ntks_Replacement(comp_ntks_cp, leaf_mapping, no_lfs, comp_ind); 
    return Extend_Final_Ntks(final_ntks_list,  comp_ntks_cpcp, node_ind);
  } /* Combine_Solutions */

void Generate_Permut(char *zh, int pos,  int new_pos, char *current_zh){
      strcpy(current_zh, zh); 
      current_zh[new_pos]=zh[pos]; 
      current_zh[pos]=zh[new_pos];
}

int Min(int a, int b){
   if (a <b) return a; else return b;
}


void Compute_Incidence_July30(int permute[], int inci[],  int  pa[], int lfs_no){
    int i,size;
    int taxa, parent;

    
    size=2*lfs_no+1;
    for (i=1; i<size; i++)  inci[i]=-1;
    inci[0]=permute[1];
    
    for (i=1; i<=lfs_no; i++) {
      taxa=permute[i];
      while (pa[taxa]!=-1) {
        parent=pa[taxa];
        if (inci[parent]==-1) { 
           inci[parent]=permute[i]; taxa=parent;
        } else { inci[parent]=permute[i]; break; }
      } /* while */
    } /* taxa */
} /* Compute_Incidence_July30 */


int ConverteIndex_Aug21(int denom[], int new_pos[], int dim){
    int result;
    int i;

    result=new_pos[0];
    for (i=1; i<dim; i++) { result=result+new_pos[i]*denom[i-1]; }
    return result;
} 


void Int_to_Table_Pos_Aug21(int num,  int dim,  int denom[], int pos[]){
     int r, q;
     int i;

     q=num;
     for (i=1; i<dim; i++) {
       pos[dim-i]=(int)(q/denom[dim-1-i]); q=q%denom[dim-1-i]; 
     } 
     pos[0]=(int)q; 
}


int Inverse(int a, int m,  int pi[]){
   int i;

   for (i=1; i<=m; i++) { if (pi[i]==a) return i; }
   return -1;
}


void Compute_Pos(int new_pos[], int pos[], int appear_index[], int dim){
    int i, num;

    for (i=0; i<dim; i++) new_pos[i]=pos[i];
    num=appear_index[0];
    for (i=1; i<=num; i++) new_pos[appear_index[i]]=pos[appear_index[i]]-1;
}


int Compute_Align_Aug29(int l_pos, int r_pos, int lfs_no, int lens[], 
   int *reduced_seqs[], int **res){

    char *score;
   int *taxa_appear[NO_LEAVES];
   int  j, k, taxa, dim;
   int new_pos[NO_LEAVES], pos[NO_LEAVES];
   int no;
   int s, s1;
   int  i, L, num, cell, opt_cell, opt_taxa;
   int seq_len;
   struct inform {
     char taxa;   /* taxa at the current position */
     int   prvs_cell; /* back track purpose */
   } *cond;
   int *denom=NULL;
   /* int mm; */


  dim=r_pos-l_pos+1;
  denom=(int *)malloc(sizeof(int)*dim);
  for (j=0; j<dim; j++) {
     if (j==0) denom[j]=(lens[j+l_pos]+1);
     else denom[j]=denom[j-1]*(lens[j+l_pos]+1);
   }

  score=(char *)malloc(denom[dim-1]+1);
  cond=(struct inform *)malloc((denom[dim-1]+1)*sizeof(struct inform));
  for (j=0; j<=lfs_no; j++)
     taxa_appear[j]=(int *)malloc(sizeof(int)*(dim+1));

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
          s1=(int)(score[cell]-'0');
         if ( s>s1  ) { s=s1;  opt_cell=cell; opt_taxa=i; }
       }
     }/* i loop */
     score[num]='0'+(s+1);
     cond[num].taxa='0'+opt_taxa; cond[num].prvs_cell=opt_cell;

  }/* num loop */
  
  num=L; seq_len=0;
  *res=(int *)malloc(sizeof(int)*(2+(score[L]-'0')));
  while (num !=0) {
      (*res)[seq_len]=(int)(cond[num].taxa-'0');
      num= cond[num].prvs_cell;
      seq_len++;
  }

  free(score); score=NULL;  free(cond); cond=NULL;
  for (j=0; j<=lfs_no; j++)
      { free(taxa_appear[j]); taxa_appear[j]=NULL;}
  return seq_len;
} /* Comput_Align_Augh29 */

/* this version only find the len of supersequecen, not itself */
int Compute_Score_HighDim(int dim,  int tip, int lfs_no, 
  int *incidence[][NO_LEAVES+1], int non0_trs1[],  int **align_res){
   char *score;
   int *taxa_appear[NO_LEAVES];
   int new_pos[NO_LEAVES];
   int pos[NO_LEAVES];
   int  no, i,  j, k, taxa;
   int  num, cell, L, S, d;
   int len, seq_len, l_pos, r_pos;
   int *res=NULL;
   int *reduced_seqs[NO_TREES];
   int lens[NO_TREES];
   int *denom, max_cells;
   /* MAX INT=2,147,483,647 */


  for (i=0; i<dim; i++){ 
       lens[i]=incidence[non0_trs1[i]][tip][0];
       reduced_seqs[i]=(int *)malloc((lens[i]+1)*sizeof(int));
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
         if (align_res!=NULL && *align_res!=NULL) free(*align_res);
         *align_res=(int *)malloc(sizeof(int)*len);
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
       reduced_seqs[r_pos]=(int *)malloc((len+1)*sizeof(int));
       reduced_seqs[r_pos][0]=len;
       for (j=1; j<=len; j++) { 
            reduced_seqs[r_pos][j]=(res)[len-j];
       }
        S=len*(lens[i]+1); l_pos=r_pos;
     } else  { 
        S=S*(lens[i]+1);
     }
   } /* for */

   for (j=0; j<dim; j++){ if (reduced_seqs[j]!=NULL) free(reduced_seqs[j]); }


} /* computer score */

 int  Compute_Score(int dim, int lfs_no,  int lens[],
  int *reduced_seqs[],  int **seq){
   char *score;
   int *taxa_appear[NO_LEAVES];
   int  j, k, taxa;
   int new_pos[NO_LEAVES], pos[NO_LEAVES];
   int no;
   int s, s1;
   int  i, L, num, cell, opt_cell, opt_taxa;
   int seq_len;
   struct inform {
     char taxa;   
     int   prvs_cell; 
   } *cond;
   int *denom=NULL;
   
  denom=(int *)malloc(sizeof(int)*(dim+1));
  for (j=0; j<dim; j++) {
     if (j==0) denom[j]=(lens[j]+1);
     else denom[j]=denom[j-1]*(lens[j]+1);
   }

  score=(char *)malloc((denom[dim-1]+1)*sizeof(char));
  cond=(struct inform *)malloc((denom[dim-1]+1)*sizeof(struct inform));
  for (j=0; j<NO_LEAVES; j++)
     taxa_appear[j]=(int *)malloc(sizeof(int)*(dim+1));

  score[0]='0';
  L=denom[dim-1]-1;
  for (num=1; num<=L; num++) {
      Int_to_Table_Pos_Aug21(num, dim, denom,  pos);

       for (j=1; j<=lfs_no; j++) { taxa_appear[j][0]=0; }
       for (j=0; j<dim; j++){
        if (pos[j]>0){  
          taxa=reduced_seqs[j][pos[j]];
          k=taxa_appear[taxa][0];
          taxa_appear[taxa][0]= k+1;
          taxa_appear[taxa][k+1]=j;
        }
      } 

     s=dim*lfs_no+1;
     for (i=1; i<=lfs_no; i++) {
       if (taxa_appear[i][0] >0) {
         Compute_Pos(new_pos, pos, taxa_appear[i], dim);
         cell=ConverteIndex_Aug21(denom, new_pos, dim);
          s1=(int)(score[cell]-'0');
         if ( s>s1  ) { s=s1;  opt_cell=cell; opt_taxa=i; }
       }
     }

     score[num]='0'+(s+1); 
     cond[num].taxa='0'+opt_taxa; cond[num].prvs_cell=opt_cell;

  }


  num=L; seq_len=0;
  if (seq!=NULL) {
     if (*seq!=NULL) { free(*seq); *seq=NULL; }
     *seq=(int *)malloc(sizeof(int)*(1+score[L]-'0'));
  }
  while (num !=0) {
      (*seq)[seq_len]=(cond[num].taxa-'0'); 
      num= cond[num].prvs_cell; 
   seq_len++;
  }

  free(score); score=NULL;  free(cond); cond=NULL;
  for (j=0; j<NO_LEAVES; j++) 
      { free(taxa_appear[j]); taxa_appear[j]=NULL;}
  return seq_len;
} 
 /* computer score */


/* two sequence alignment */
int Align_TwoIntSeq(int seq1[], int seq2[]){
   int i, j, k;
   /*
   int score[NO_NODES][NO_NODES];
   */
   int len1, len2, len;
   int *score[NO_NODES];

   len1=seq1[0]; len2=seq2[0];
   for (i=0; i<=len1; i++) {
       score[i]=(int *)malloc(sizeof(int)*(len2+1));
   }

   for (i=0; i<=len1; i++) {
     for (j=0; j<=len2; j++) {
        if (i==0) score[i][j]=j;
        else if (j==0) score[i][j]=i;
        else if (seq1[i]==seq2[j]) score[i][j]=1+score[i-1][j-1];
        else score[i][j]=1+Min(score[i][j-1], score[i-1][j]);
     }
   } /* i loop */
 
   k=score[len1][len2];
   for (i=0; i<=len1; i++) free(score[i]);
   return k;
}/* Align_TwoIntSeq( */

/* we assume the fist list is inter */
int Is_Subsequence(int list1[], int list2[]){
  int i, s1, s2;
  int k;

  s1=list1[0]; s2=list2[0];
  k=1;
  for (i=1; i<=s2; i++) {
    if (list2[i]==list1[k]) k++;
    if (k==s1+1) return 1;
  }
  
  return 0;

}

int Equal_Containment_Aug21(int list1[], int list2[]){
   int i, j, k;

  if (list1[0]<list2[0]) { 
        return Is_Subsequence(list1, list2);
  } else {
        return Is_Subsequence(list2, list1);
  }

}

/* return o if it is not contained in any sequence after that */
int Contain(int *incidence[][NO_LEAVES+1], int non0_trs[],
   int  id,  int  dim, int taxa){

  int  i;

  for (i=id; i<dim; i++) {
    if (Equal_Containment_Aug21(incidence[non0_trs[id-1]][taxa],
          incidence[non0_trs[i]][taxa])==1) 
     return 1;
  } 
 return 0;

} 


void Convert_Into_Array(int *inci_str[], int taxa_inci[], int lf_no, 
 int p[]) {

    int i, k; 
    int pa, j;
    int  str[NO_LEAVES];

    /* printf("       inside Convert\n"); */
    inci_str[0]=NULL;
    for (i=1; i<=lf_no; i++) {
      pa=p[i]; j=0;
      if (pa!=-1) {
        while (taxa_inci[pa]!=i) {
           str[j]=taxa_inci[pa]; 
           j++; pa=p[pa];
        }
      }
      inci_str[i][0]=j;
      for (k=1; k<=j; k++) inci_str[i][k]=str[k-1];
    }
}

int Is_Contain_Str(char *str,  char *sub_str) {
    int i, s1, s2;
  int k;

  s2=strlen(str);
  s1=strlen(sub_str);
  k=0;
  for (i=0; i<s2; i++) { 
      if (k<s1 && sub_str[k]==str[i]) k++; 
      if (k==s1) return 1; 
  }

  return 0;

}

int  Equality_Contain_str_200(char seqs_str[][NO_NODES], int k,  char *str){ 
     int i;

     for (i=0; i<=k; i++) {
        if (Is_Contain_Str(seqs_str[i], str) ==1) return 1;
     }
     return 0;
}


int Subsequence_Check_100(char seqs_str[][NO_NODES], int no,  char *str){

     int i;
     for (i=0; i<no; i++) {
        /* if (strcmp(seqs_str[i], str)==0  ) return 1; */
        if (Is_Contain_Str(str, seqs_str[i]) ==0) return 0;
     }
     return 1;
}


int Subsequence_Check(char *seqs_str[], int no,  char *str){

     int i;
     for (i=0; i<no; i++) {
        /* if (strcmp(seqs_str[i], str)==0  ) return 1; */
        if (Is_Contain_Str(str, seqs_str[i]) ==0) return 0;
     }
     return 1;
}

int Reduction_100(char seqs_str[][NO_NODES], int no, char *align_result){
   int len;
   int i, j;
   char x;
   char *str;

   len=strlen(align_result); 
str=(char *)malloc(strlen(align_result)+1);
   strcpy(str, align_result);

   for (i=0; i<len; i++) {
    str[i]='.';
    if (Subsequence_Check_100(seqs_str, no, str)==0) str[i]=align_result[i];
   }


   i=0; j=0;
   while (j<len) {
     while (str[j]=='.') j++;
     if (j<len) {
      align_result[i]=str[j];
      i++; j++;
     }
   }

   align_result[i]='\0';
   free(str);
   return i;
}


int Reduction(char *seqs_str[], int no, char *align_result){
   int len;
   int i, j;
   char x;
   char *str;

   len=strlen(align_result); 
   str=(char *)malloc(strlen(align_result)+1);
   strcpy(str, align_result);

   for (i=0; i<len; i++) {
    str[i]='.';
    if (Subsequence_Check(seqs_str, no, str)==0) str[i]=align_result[i];
   }


   i=0; j=0;
   while (j<len) { 
     while (str[j]=='.') j++;
     if (j<len) {
      align_result[i]=str[j];
      i++; j++;
     }
   }

   align_result[i]='\0';
   free(str);
   return i;
}


int Align2_July8(int no_trees, int taxa, 
  int *incidence[][NO_LEAVES+1], int lfs_no, int **align_res){

  int i, j, m;
  int len[NO_TREES], non0_trs[NO_TREES], pos[NO_TREES];
  int non0_trs1[NO_TREES];
  int min;
  int dim,  dim1, index;
  int seq_len, seq_len2;

  char *seq_str1[NO_TREES];
  char *align_res_str;
 /*
  int *reduced_seqs[NO_TREES]={NULL};
  */

  dim=0; 
  for (i=0; i<no_trees; i++) {
     len[i]=(-1)*incidence[i][taxa][0];
     if (incidence[i][taxa][0]>0) { non0_trs[dim]=i; dim++; }
  } /* non0_trs contains all trees where the first components >0)*/


  if (dim==0) { *align_res=NULL;  return 0; }
  else if (dim==1) {
    seq_len=incidence[non0_trs[0]][taxa][0];
    *align_res=(int *)malloc(sizeof(int)*(seq_len)); 
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
      *align_res=(int *)malloc(sizeof(int)*(seq_len+1)); 
      for (j=1; j<=seq_len; j++) 
        (*align_res)[j-1]=incidence[non0_trs1[0]][taxa][seq_len-j+1];  
      return seq_len;
  }

   dim1=dim;
   for (i=0; i<dim; i++) {
     seq_len=len[non0_trs1[i]]*(-1);
     seq_str1[i]=(char *)malloc(seq_len+1);
     for (j=1; j<=seq_len; j++) 
        seq_str1[i][j-1]=incidence[non0_trs1[i]][taxa][seq_len-j+1];
     seq_str1[i][seq_len]='\0';
   }

  seq_len=Compute_Score_HighDim(dim,taxa,lfs_no,incidence, non0_trs1,align_res);

   align_res_str=(char *)malloc(seq_len+1);
   for (i=0; i<seq_len; i++) align_res_str[i]='0'+((*align_res)[i]-1);
    align_res_str[seq_len]='\0';
     
    seq_len2=Reduction(seq_str1, dim1, align_res_str);
    if (seq_len2<seq_len) {
      /* a potential problem */
       for (j=0; j<seq_len2; j++) 
         (*align_res)[j]=1+(align_res_str[j]-'0');
    }

    for (i=0; i<dim1; i++) { if (seq_str1[i]!=NULL) free(seq_str1[i]);}
  /*
  printf("   Align2_July out\n");
  */

  return seq_len2;
  
} /* Align2_July8 */

int Pairwise_Align_July8(char *seq1, char *seq2){
   int i, j;
   int **score;
   int len1, len2, len;

   len1=strlen(seq1); len2=strlen(seq2);
   score=(int **)malloc(sizeof(char *)*(len1+1));
   for (i=0; i<=len1; i++) score[i]=(int *)malloc(sizeof(int)*(len2+1));

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
int Pairwise_Align(char *seq1, char *seq2, char seq_str[]){
   int i, j, k;
   int **score;
   int len1, len2, len;
   char *seq;
   char t;

   /* printf("inside: PW align\n"); */
   len1=strlen(seq1); len2=strlen(seq2);
   score=(int **)malloc(sizeof(char *)*(len1+1));
   for (i=0; i<=len1; i++) score[i]=(int *)malloc(sizeof(int)*(len2+1)); 

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
   free(seq);
   return len;
}


/*
int Align3_1000(int no_trees, int taxa, char *incidence[][NO_LEAVES], int lfs_no, int **seq, int checker){
*/
int Align3_1000(int no_trees, int taxa, char **incidence[], int lfs_no, int **seq, int checker){

  int i, j, k, m;
  int  min, dim,  dim1;
  int seq_len, seq_len2;
  int q_flag, coord[TWO], co[TWO];
  char seqs_str[NO_TREES][NO_NODES];
  char seqs_str1[NO_TREES][NO_NODES];
  char seqs_str3[NO_TREES][NO_NODES];
  char  *align_result=NULL;
  int *len=NULL, *array=NULL;
  int **scores=NULL;

/* Align the LTS of Taxa in the different trees, 
 * do multiple sequence alignemtn in the limited space */


  dim1=0;
  for (i=0; i<no_trees; i++) {
     k=strlen(incidence[i][taxa]);
     if (k>0) {
       strcpy(seqs_str3[dim1], incidence[i][taxa]);
       dim1++;
     } 
  } 

  if (dim1==0) {
     if (seq!=NULL && (*seq!=NULL) ) { free(*seq); *seq=NULL;} 
     return 0;
  }

  /* convert a leter into taxon index */
  if (dim1==1) {
      seq_len=strlen(seqs_str3[0]);
      if (seq!=NULL) {
         if (*seq!=NULL) {free(*seq); *seq=NULL;}
         *seq=(int *)malloc(sizeof(int)*seq_len);
         for (j=0; j<seq_len; j++) (*seq)[j]=1+(seqs_str3[0][seq_len-1-j]-'0');
      }
      /*
      free(seqs_str1[0]);  
      */
      return seq_len;
  }


  len=(int *)malloc(dim1*sizeof(int));
  array=(int *)malloc(dim1*sizeof(int));
  for (i=0; i<dim1; i++) { len[i]=strlen(seqs_str3[i]); array[i]=i; }
  MergeSort(array, 0, dim1-1, len);



  if (checker==1) {
    for (i=0; i<dim1; i++) printf("array[%d]=%d\n", i, array[i]);
  }
  
  i=0; j=1;
  strcpy(seqs_str[0], seqs_str3[array[0]]);

  /* remove identical copies */
  while (j<dim1) {
    if (Equality_Contain_str_200(seqs_str, i,  seqs_str3[array[j]])==0){
      i=i+1;
      strcpy(seqs_str[i], seqs_str3[array[j]]);
    }
    j=j+1;
  }
  dim=i+1;
  free(len);
  free(array);



  if (dim==1) {
     seq_len=strlen(seqs_str[0]);
     if (seq!=NULL) {
         if (*seq!=NULL) {free(*seq); *seq=NULL;}
         *seq=(int *)malloc(sizeof(int)*seq_len);
         for (j=0; j<seq_len; j++) (*seq)[j]=1+(seqs_str[0][seq_len-1-j]-'0');
     }
     return seq_len;
 }


  scores=(int **)malloc(sizeof(int *)*(dim+1));
  for (i=0; i<dim; i++) {
       scores[i]=(int *)malloc(sizeof(int)*(dim+1));
       strcpy(seqs_str1[i], seqs_str[i]);
  }
  dim1=dim;


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
       if (q_flag==0) { min=scores[i][j]; coord[0]=i; coord[1
]=j; q_flag=1;
       } else {
         if (min > scores[i][j]) { min=scores[i][j]; coord[0]
=i; coord[1]=j; }
       }
     }
   }

  align_result=(char *)malloc(strlen(seqs_str[coord[0]])+strlen(seqs_str[coord[1]])+1);
 seq_len=Pairwise_Align(seqs_str[coord[0]], seqs_str[coord[1]], align_result);

   strcpy(seqs_str[coord[0]], align_result);
  free(align_result);

   if (coord[1] !=dim-1) {
      strcpy(seqs_str[coord[1]], seqs_str[dim-1]);
   }
   co[0]=coord[0]; co[1]=coord[1];
   dim--;
  } /* while */

  /* finding the supersequence */
   align_result=(char *)malloc(strlen(seqs_str[0])+strlen(seqs_str[1])+1);
     seq_len=Pairwise_Align(seqs_str[0], seqs_str[1], align_result);
     seq_len2=Reduction_100(seqs_str1, dim1, align_result);

  if (seq!=NULL) {
     if (*seq!=NULL) { free(*seq); *seq=NULL;}
     *seq=(int *)malloc(sizeof(int)*strlen(align_result));
     for (j=0; j<seq_len2; j++) (*seq)[j]=1+(align_result[seq_len2-1-j]-'0');
     /* free(align_result); */
  }


  for (i=0; i<dim1; i++) {
    /*
    if (seqs_str[i]!=NULL) free(seqs_str[i]);
     */
       /*  if (seqs_str1[i]!=NULL)   free(seqs_str1[i]); */  
      free(scores[i]);
    }
  free(scores);
  return seq_len2;
} /* Align3_1000 */

  

/* in the following construction, the leaves are 1 to lfs_no;
 * the potential ret are lfs_no+2,..., 2lfs_no. tree nodes are
 *  2lfs_no+1... 
 */
/* construct tree-child ntks and return the no. of edges */
int Compute_NTK(int no_trees, int lfs_no, int ret_bound,
 int permut[], int  *inci[][NO_LEAVES+1], 
  char  **incidence[], int *gph1[2], int *no_nodes){
     int a, b, i, j, k, m, u, v;
     int seq_len, seq_len1;
     int no_edges, node;
     int ret;
     int *seq=NULL;
     int *indeg=NULL;
     int *parent=NULL;
     int *gph[TWO];
     int no_deg1_nodes;
     int **align_res=NULL;

   /*
     printf("in: comput ntk, %d \n", 3*lfs_no+2*ret_bound);
    */

    // gph[0]=(int *) malloc(sizeof(int)*(3*lfs_no+2*ret_bound));
    // gph[1]=(int *) malloc(sizeof(int)*(3*lfs_no+2*ret_bound));
     gph[0]=(int *) malloc(sizeof(int)*MAXSIZE);
     gph[1]=(int *) malloc(sizeof(int)*MAXSIZE);

    no_edges=0; node=2*lfs_no+1;

   if (align_res!=NULL) { free(align_res); align_res=NULL;}
         align_res=(int **)malloc(sizeof(int *));


    for (i=1; i<lfs_no; i++) {
        m=permut[i];
       
       /* high dimension alignment */ 
       seq_len=Align2_July8(no_trees, m, inci,  lfs_no, align_res);
       if (no_trees > 2){ 
          /* greedy apprach */
          seq_len1=Align3_1000(no_trees, m, incidence,  lfs_no, &seq, 0);
          if (seq_len1 < seq_len) { 
            seq_len=seq_len1;
            if (*align_res!=NULL) { free(*align_res); *align_res=NULL;}
            *align_res=(int *)malloc(sizeof(int)*seq_len);
            for (j=0; j<seq_len; j++)  (*align_res)[j]=seq[j];
          }
        }

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

    if (align_res!=NULL) { 
     if (*align_res!=NULL) free(*align_res);
     free(align_res); align_res=NULL;
    }

    
    parent=(int *)malloc(sizeof(int)*(2*lfs_no+2));
    indeg=(int *)malloc(sizeof(int)*(2*lfs_no+2));
    /* remove degree-two nodes */ 
    a=lfs_no+1; b=2*lfs_no+1;
    for (i=a; i<b; i++) { indeg[i]=0; } 
    for (j=0; j<no_edges; j++) {
       m=gph[1][j]; 
       if (m>=a && m<b) {
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

    
    free(indeg); free(parent);
    return k;
}


/* return ret number */
int Search_minNTKS_Round1(int no_trees, int lfs_no, int permut[], 
char **incidence[], int pos, int checker){
     int i;
     int  seq_len, no_rets, no_ret_edges;
     int *tmp=NULL;
     

    no_ret_edges=0; 
    
    for (i=1; i<=pos; i++) {
    seq_len=Align3_1000(no_trees, permut[i], incidence,  lfs_no, &tmp, checker);
       no_ret_edges=no_ret_edges+seq_len;
    }

    no_rets=no_ret_edges-pos;
    if (tmp!=NULL) {free(tmp); tmp=NULL;}
    return no_rets;
    
} /* search min_ntks */

/* new version we use the exact align score */
int Search_minNTKS_Round2(int no_trees, int lfs_no, int permut[], 
int  *inci[][NO_LEAVES+1], char **incidence[], int pos){
     int i, j;
     int  seq_len, seq_len1;
     int no_rets, no_ret_edges;
     int **align_res=NULL;
     int *seq=NULL;

    
    no_ret_edges=0; 
    align_res=(int **)malloc(sizeof(int *));
    for (i=1; i<pos; i++) {
       seq_len=Align2_July8(no_trees, permut[i], inci,  lfs_no, align_res);
       if (no_trees >2 ){
          seq_len1=Align3_1000(no_trees, permut[i], incidence,  lfs_no, &seq, 0);
          if (seq_len1 < seq_len) seq_len=seq_len1;
       }
       
       if (*align_res!=NULL) {free(*align_res); *align_res=NULL;}
       no_ret_edges=no_ret_edges+seq_len;
    }

    no_rets=no_ret_edges-(pos-1);
    free(align_res); 
    return no_rets;
    
} /* search min_ntks */



int Is_Tree_Node(int in_deg){
     if (in_deg==1) return 1; else return 0;
}


int Is_in_Array(int a,  int k2, int array[]){
   int i;
   for (i=1; i<=k2; i++) {
     if (a==array[i]) return 1;
   }
   return 0;
}

/* we assume parents are sorted from smallest to the larger */
int Is_Equal(int array1[], int array2[]){
   int i, j, k1, k2;

   k1=array1[0]; k2=array2[0];
   if (k1!=k2) return 0;
   else {
     for (i=1; i<=k1; i++) {
        if (Is_in_Array(array1[i], k2, array2)==0) return 0;
     }
     return 1;
   }
}

int Is_Isomorphism_Dec1(int sorted1[], int *p1[],
  int *p2[], int no_nodes, int no_rets, int no_lfs){

   int i, j, k; 
   int *map, *parents;
   int xxx, node, image, node2;
   int parent1, parent2;



   map=(int *)malloc(sizeof(int)*(no_nodes+1));
   parents=(int *)malloc(sizeof(int)*(no_nodes+1));

   map[0]=0;
   for (i=1; i<=no_lfs; i++) { map[i]=i;}
   for (i=no_lfs+1; i<=no_nodes; i++) { map[i]=0;}

   for (i=no_nodes-1; i>0; i--) {
       node=sorted1[i];
      if ( Is_Tree_Node(p1[node][0])==1) { /* tree node */
         image=map[node];
         parent1=p1[node][1];
         parent2=p2[image][1];

         if (map[parent1]==0){
           if ( Is_Tree_Node(p1[parent1][0])==1){ /* parent is trnode*/
                if  ( Is_Tree_Node(p2[parent2][0])==1) {
                    map[parent1]=parent2;
                 }  else return 0;
           } else {
               if  (Is_Tree_Node(p2[parent2][0])==1) return 0;
               else  { map[parent1]=parent2; }
           }
        } else {
          if (map[parent1]!=parent2) return 0;
        }
     } /* if */
   } /* for */
 

   for (i=no_nodes-1; i>0; i--) {
       node=sorted1[i];
       node2=map[node];
       if (Is_Tree_Node(p1[node][0])==0){
          parents[0]=p1[node][0];
          for (j=1; j<=p1[node][0]; j++)  parents[j]=map[p1[node][j]];
          if (Is_Equal(parents, p2[node2])==0) return 0;
       }
   }

   free(map); free(parents);
   
   return 1;
} /* Dec1 */


int  Check_Iso_v2(int *gph[], g_set *ntks, int edge_no, 
   int node_no, int rets_no,  int lfs_no, int *pa[]){
   /*
   int  map[NO_NODES];
   */
   g_set *ntk_ptr;
   int i, j,k, ans;
    int **pa2; 
   int *sorted_nodes=NULL;
   int  counter;



      pa2=(int **)malloc(sizeof(int *)*(node_no+1));
      sorted_nodes=(int *)malloc(sizeof(int)*(node_no+1));

      top_sort_Aug24(gph,  node_no, edge_no, sorted_nodes);

     ntk_ptr=ntks;
     counter=0;
     while (ntk_ptr!=NULL) {
       counter +=1;
       if (node_no==ntk_ptr->no_nodes){
         pa2[0]=(int *)malloc(sizeof(int));
         pa2[0][0]=(ntk_ptr->parent)[0][0];
        for (i=1; i<node_no+1; i++) {
          k=(ntk_ptr->parent)[i][0];
          pa2[i]=(int *)malloc((k+1)*sizeof(int));
          for (j=0; j<=k; j++) {
              pa2[i][j]=(ntk_ptr->parent)[i][j];  
          }
        }
        ans=Is_Isomorphism_Dec1(sorted_nodes, pa,  pa2, node_no, rets_no, lfs_no);
        for (i=0; i<node_no+1; i++) { free(pa2[i]); }
        if (ans==1) { 
           free(sorted_nodes); 
           free(pa2);
           return 1; 
        }
      } else { }
      ntk_ptr=ntk_ptr->next;
   } /* while */

   free(sorted_nodes);
   free(pa2);
   return 0;
}

g_set *Insert_v2(g_set *ntk_list, int *gph[TWO], int edge_no, int node_no,
  int ret_no, int lfs_no, char *zh){
   g_set *ptr=NULL;
   int i,k,j, end;
   int ans;
   int *pa[NO_NODES];
   
   Compute_parents(gph, edge_no, node_no,  pa);

    if (ptr!=NULL) {free(ptr); ptr=NULL;}
    if (ntk_list==NULL) {
       ptr=(g_set *)malloc(sizeof(g_set));
       ptr->no_nodes=node_no;
       ptr->no_edges=edge_no;
       ptr->left_ends=(int *)malloc(edge_no*sizeof(int));
       ptr->right_ends=(int *)malloc(edge_no*sizeof(int));
       for (i=0; i<edge_no; i++) {
         (ptr->left_ends)[i]=gph[0][i+1];
         (ptr->right_ends)[i]=gph[1][i+1];
       }
      /* (ptr->parent)=pa; */
       (ptr->parent)[0]=(int *)malloc(sizeof(int));
       (ptr->parent)[0][0]=-1;
       for (i=1; i<node_no+1; i++) { 
           k=pa[i][0];
           (ptr->parent)[i]=(int *)malloc((k+1)*sizeof(int));
           for (j=0; j<=k; j++) (ptr->parent)[i][j]=pa[i][j];
       }
       ptr->next=NULL;
       for (i=0; i<=node_no+1; i++) { if (pa[i]!=NULL) free(pa[i]);}
       return ptr;
     } else {
      ans=Check_Iso_v2(gph, ntk_list, edge_no, node_no, ret_no, lfs_no, pa);

      if (ans==0) {
           /* printf("permutation: %s\n", zh); */
           ptr=(g_set *)malloc(sizeof(g_set));
           ptr->no_nodes=node_no;
           ptr->no_edges=edge_no;
           ptr->left_ends=(int *)malloc(edge_no*sizeof(int));
           ptr->right_ends=(int *)malloc(edge_no*sizeof(int));
           for (i=0; i<edge_no; i++){ 
                 (ptr->left_ends)[i]=gph[0][i+1];
                 (ptr->right_ends)[i]=gph[1][i+1];
           }
           (ptr->parent)[0]=(int *)malloc(sizeof(int));
           (ptr->parent)[0][0]=-1;
           for (i=1; i<node_no+1; i++) { 
               k=pa[i][0];
              (ptr->parent)[i]=(int *)malloc((k+1)*sizeof(int));
              for (j=0; j<=k; j++) (ptr->parent)[i][j]=pa[i][j];
            }
            ptr->next=ntk_list;
    for (i=0; i<=node_no+1; i++) { if (pa[i]!=NULL) free(pa[i]); }  
            return ptr;
      } else {  
    for (i=0; i<=node_no+1; i++) { if (pa[i]!=NULL) free(pa[i]); }  
       return ntk_list; }
   }/* if */

   
}

int Comput_Cutoff(int ntk_counts[], int min, int max, int total){
    int i, j;
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
  int cut_off,  int min_ret_no, int ctt){

    struct node *Qptr, *f_ptr, *e_ptr, *bad;

    
     /* srand( time(0)); */     
     srand(123);
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

void Convert_Str_Format_July30(char *inci_str[],int taxa_inci[],int lf_no,int p[]){
    int i; 
    int pa, j;
    char *str;


    str=(char *)malloc(lf_no+1);
 
    for (i=1; i<=lf_no; i++) {
      pa=p[i]; j=0;
      /* while (taxa_inci[pa]!=i) { */ 
      if (pa!=-1) {
        while ( 1 ) { 
          if (taxa_inci[pa]==i) break;
          str[j]='0'+(taxa_inci[pa]-1);
          j++; pa=p[pa];
        }
      }
      
      str[j]='\0';
   
      strcpy(inci_str[i], str);
    }
   free(str);
} /*  Convert_Str_Format_July30( */


int  Find_Id(char *taxa_names[], int index,  char *p1){
   int i, j;

   for (i=0; i<index; i++) {
     if (strcmp(taxa_names[i], p1)==0)  return i;
   }
   return -1;
}



void  Replace_Taxa_with_ID_July20(char* str, char* taxa_names[],  int index) {
    char* p = NULL;
    int i,id,k;
    int len = 0;
    char new_str[TREE_INPUT_STR]={0};
    char id_str[TAXA_NAME_STR]={0};
    char p1[TAXA_NAME_STR]={0}, p2[TAXA_NAME_STR]={0};
    int low, high;

    len=strlen(str);
    high=0; id=-1; low=0;
    k=0;
   for (high = 0; high <len; high++) {
        if (str[high] == ',' || str[high] == ')') {
            if (high != low) {
                memcpy(p1, &str[low], high - low);
                p1[high - low]='\0';
                id=Find_Id(taxa_names, index, p1);
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
        }
   }
   strcpy(str, new_str);
}


/* number the internal nodes */
void Assign_InternalIDs(char* c, int taxa_no) {
    int i, len = 0;
    int n = 0, amount;
    char a[TREE_INPUT_STR];

    amount=taxa_no+2;
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
int  Data_search(char* nw_str, int i, int flag, int l, 
   char* str, int tree[NO_NODES][2], int current_edge) {

    int j = 0, k;
    char taxa_id2[TAXA_NAME_STR] = { 0 }, taxa_id3[TAXA_NAME_STR] = { 0 };
    char a[TAXA_NAME_STR] = { 0 };

    if (flag) {
        while (nw_str[i] >= 48 && nw_str[i] <= 57 && i < l)
        {
          sprintf(taxa_id2, "%d", nw_str[i] - 48); 
          strcat(taxa_id3, taxa_id2);
          i++;
        }
        tree[current_edge][0] = (int) atoi(taxa_id3);
        /* (c1 -> flag3)++; */
    } else {
        while (str[i] >= 48 && str[i] <= 57) { a[j++] = str[i--]; }
        for (k = strlen(a) - 1; k >= 0; k--) {
            sprintf(taxa_id2, "%d", a[k] - 48); strcat(taxa_id3, taxa_id2);
        }
        tree[current_edge - 1][1] = (int) atoi(taxa_id3);
        return i + 1;
    }
}


/*   edge_no=Compute_tree(newick_str, tree, taxa_names, index); */
/*Connect parent and child nodes to generate tree*/
int  Compute_tree(char* nw_str, int tree[NO_NODES][2]) {
    char a[TAXA_NAME_STR] = { 0 };
    int i, j, current_edge;

   /* printf("Comput_tree: %d   %s\n", index, nw_str); */

    j=0;
    current_edge=0;
    for (i = 0; i < strlen(nw_str); i++) {
        a[j++] = nw_str[i];
        if (nw_str[i] == ')') {
            while (a[--j] != '(') {
              Data_search(nw_str, i + 1, 1, strlen(nw_str), a, tree, current_edge);
                current_edge++;
              j = Data_search(nw_str, j - 1, 0, strlen(nw_str), a, tree, current_edge);
        /*
        printf("current_edge(%d), left(%d %d)\n", current_edge, tree[current_edge-1][0], tree[current_edge-1][1]);
         */
            }
       }
    }
    return  current_edge;
}


int Is_New_Names(char *taxon_name, char *names[], int old_index){
      int k;

      /*
      printf("inside(Is_new: old_index(%d); taxon_name(%s)\n", old_index, taxon_name); 
       */
      if (old_index==0) return 0;
      for (k=0; k<old_index; k++) {
        /*
        printf("taxon_name(%s); ", taxon_name); 
        printf(" names_array(%s)\n",  names[k]); 
        */
        if (strcmp(taxon_name, names[k])==0) { return 1; }
      }
      return 0;

}


int Extract_Taxa(char *newick_str, char *names[], int *current_index){
    int i;
    char taxa_names[NO_LEAVES][TAXA_NAME_STR];
    int  index, low, high, len;
    int old_index;

    len=strlen(newick_str);
    if (newick_str[len - 1] != ')') { newick_str[len - 1] = '\0'; }
    high=0; index=0;  old_index=*current_index; low=0;
    for (high = 0; high <len; high++) {
        if (newick_str[high] == ',' || newick_str[high] == '(' ||
           newick_str[high] == ')') {
            if (high != low) {
                memcpy(taxa_names[index], &newick_str[low], high - low);
                /* printf("index(%d), new_taxa(%s)\n", index, taxa_names[index]); */
                if (Is_New_Names(taxa_names[index], names,  old_index)==0){
                  /*  printf(" here\n"); */
                 names[old_index]=(char *)malloc(strlen(taxa_names[index])+1);
                 strcpy(names[old_index], taxa_names[index]);
                 old_index++;
                }
                index++;
            }
            low = high+1;
        }
       /*  high++; */
   }
   *current_index=old_index;
   return index;
} /* Extract_Taxa */


int  Convert_to_Edges(char *newick_str, char *taxa_names[],
    int tree[][2], int *taxa_no, int lfs_no) {

    int  edge_no;

   Extract_Taxa(newick_str, taxa_names, taxa_no);

   /* printf("Old_tring: %s\n", newick_str); */
   Replace_Taxa_with_ID_July20(newick_str, taxa_names, *taxa_no);
   /* printf("new_tring: %s-- %d\n", newick_str, *taxa_no); */
   Assign_InternalIDs(newick_str, lfs_no);
  /*  printf("after add_intenral new_tring: %s\n\n", newick_str); */
   edge_no=Compute_tree(newick_str, tree);

   return edge_no;
} /* Convert_to_Edges */



void Add_Rootedge(int edge_no, int  tree[][2], int taxa_no){
   int node_no;
   int *in_deg, *out_deg;
   int i, root;

  /*  printf("inside_AddRoot---%d--\n", edge_no);*/
   in_deg=(int *)malloc((2*taxa_no)*sizeof(int));
   out_deg=(int *)malloc((2*taxa_no)*sizeof(int));

   for (i=0; i<= 2*taxa_no; i++) {
     in_deg[i]=0; out_deg[i]=0;
   }

   for (i=0; i< edge_no; i++) {
     in_deg[tree[i][1]] +=1; out_deg[tree[i][0]] +=1;
   }

   for (i=0; i<= 2*taxa_no; i++) {
     if (in_deg[i]==0 && out_deg[i]==2) {root=i; break;}
   }

   tree[edge_no][0]=tree[0][0]; tree[edge_no][1]=tree[0][1];
   tree[edge_no+1][0]=tree[1][0]; tree[edge_no+1][1]=tree[1][1];
   tree[1][0]=0; tree[1][1]=root;
} /* Add_Rootedge */  

 int Look_index(int node, int sorted_nodes[], int no_nodes){
   int i;
   for (i=0; i<no_nodes; i++){
     if (sorted_nodes[i]==node) return i;
   } 
 } 


void   Print_Ntks_500(ntks_struct *new_ntks, int rets, char *taxa_names[], int no_lfs, FILE *out_file){
       int graph[MAXSIZE][2];
       int i,  k,j;
       ntks_struct *ptr;
       int no_nodes, no_edges;
       int no_ntks=0;
       int *sorted_nodes=NULL;
       int ans; 



   /*  printf("---inside --Print_Ntks--\n"); */
   if (out_file ==NULL) 
     printf(" ===it is nULL==\n");

       ptr=new_ntks; 

       while (ptr!=NULL) {
         no_ntks +=1;
         no_edges=ptr->no_edges;
         
         for (j=0; j<no_edges; j++) {
               graph[j+1][0]=ptr->lt_ends[j]; graph[j+1][1]=ptr->rt_ends[j];
               /*
               printf("%d %d\n", graph[j+1][0], graph[j+1][1]);
               if (j==no_edges) printf("\n End\n");
                */
         }
            
         no_nodes=Simplify_Ntks_Aug22(graph, no_edges);
             sorted_nodes=(int *)malloc(sizeof(int)*(no_nodes+1));
     top_sort_Dec25(graph, no_nodes, no_edges, no_lfs, sorted_nodes);
     /*
     fprintf(out_file, "=============// %d, %d; \n", no_lfs, no_nodes);
      */
    

            /* fprintf(out_file, "new_form\n"); */
     for (j=1; j<=no_edges; j++) {
          if (graph[j][0]>no_lfs) 
           graph[j][0]=no_nodes-1-Look_index(graph[j][0], 
             sorted_nodes, no_nodes-no_lfs-1);
          else if (graph[j][0]>0 && graph[j][0]<= no_lfs) {
          }
 
         if (graph[j][1]>no_lfs) 
            graph[j][1]=no_nodes-1-Look_index(graph[j][1], 
                          sorted_nodes, no_nodes);
         else if (graph[j][1]<= no_lfs) {
         }
             
         if ( graph[j][1] <=no_lfs)
         fprintf(out_file, "%d %s\n", graph[j][0], taxa_names[graph[j][1]-1]);
         else 
         fprintf(out_file, "%d %d\n", graph[j][0], graph[j][1]);
         /* printf("%d %d\n", graph[j][0], graph[j][1]); */
     }

          ptr=ptr->next;
     } /* while */

      for (i=1; i<=no_lfs; i++) 
        fprintf(out_file, "%i ==> %s\n", i, taxa_names[i-1]);

   printf("\n\n    %d ntks  with %d rets are found\n", no_ntks, rets);

} /* print_nks_1 */


void Print_input_trees(int t, int no_edges, int tree[][TWO]){
   int i, j, k;

       printf(" %d-th tree with %d edges:\n", t, no_edges);
     for (j=1; j<=no_edges; j++) {
       printf("%d  %d\n", tree[j][0], tree[j][1]);
     }
     printf("\n");
}

int Not_All_Empty(int no_trees, int I_trees[][NO_NODES][TWO]){
    int i;
   for (i=0; i<no_trees; i++) {
     if (I_trees[i][0][0]>1) return 1;
   }
   return 0;
}


void Move_Largest_Tree(int no_trees,  int I_trees[][NO_NODES][TWO]){
   int i, index, size;
   int tree[NO_NODES][TWO];
   size=0; index=-1;
   for (i=0; i<no_trees; i++) {
     if (I_trees[i][0][0]>size) {index=i; size=I_trees[i][0][0];}
   }
   for (i=1; i<=size; i++) { 
     tree[i][0]=I_trees[index][i][0]; tree[i][1]=I_trees[index][i][1]; 
   }
   I_trees[index][0][0]=I_trees[0][0][0];
   for (i=1; i<=I_trees[0][0][0]; i++) {
      I_trees[index][i][0]=I_trees[0][i][0];
      I_trees[index][i][1]=I_trees[0][i][1];
   }
   for (i=1; i<=size; i++) {
     I_trees[0][i][0]=tree[i][0]; I_trees[0][i][1]=tree[i][1];
   }
   I_trees[0][0][0]=size;
}

void main(int argc, char *argv[]){

FILE *input_file, *ntk_input_file,  *tr_input_file, *ntk_file, *cluster_file; 
char *node_clusters[NO_TREES][NO_NODES];
int node1, node2, size_xxx;

int i, j, k, k1;
int no_rets, no_lfs, no_ntk_edges, no_ntk_nodes; /* ntkwork parameters */
int no_trees; 
int subtrees[NO_TREES][NO_NODES][TWO];

int I_trees[NO_TREES][NO_NODES][TWO]; /* keep ith tree in I_trees[i] */
int I_lf_no, I_edge_no, I_node_no;
int I_tree_ch[NO_TREES][NO_NODES][MAX_OUT];

int subtr_edge_no, subtr_lf_no;
int leaf_mapping[NO_TREES][NO_LEAVES+1];
int top_leafmap[NO_LEAVES+1];
char *top_clusters[NO_LEAVES+1];
int ComponentGraph[NO_NODES][NO_LEAVES]; 
int node_ind, current_comp;
int current_size;
int t_ch[NO_TREES][NO_NODES][MAX_OUT], *t_p[NO_TREES];
/*
int tr_sorted[NO_TREES][NO_NODES];
 */
int sorted_tnodes[NO_TREES][NO_NODES];
int no_components;

struct node *Q_front, *Q_end, *Qptr;
struct node  *ptr;
struct node **Second_Q[TWO];
int q_flag, min_ret_no, new_min_ret_no, current_ret_no; 
int min_node_no, ret_val;
int permut[NO_LEAVES+1];
char zh[NO_LEAVES+1], new_zh[NO_LEAVES+1];
int *gph[TWO], no_gph_edges, no_gph_nodes;
int *taxa_inci[NO_TREES];
char **taxa_incidence_str[NO_TREES];
 int *taxa_inci_int[NO_TREES][NO_LEAVES+1];
      /* indicate in the tree decomposition, whihc taxa is right below which
         and the depth of their least common ancestor */
int ctt;
int max, kkk, cut_off;
int *ntks_counts=NULL;
int  it_times, pos;
int total_rets=0;
int  xxx;
long int  count_ntks;

char input_nw_str[TREE_INPUT_STR];
char *taxa_names[NO_LEAVES];
int checker;

int mix_edge_no, taxa_no=0;

int no_lfs_top=0; 

ntks_struct  *final_ntks_list;
g_set *ntk_list;
/* INT_MAX=2,147,483,647 */


   if (argc>=4) {
      tr_input_file=fopen(argv[1], "r");
      I_lf_no=atoi(argv[2]);
   ntk_file=fopen(argv[3], "w");
      if (tr_input_file==NULL ) { printf("File openning fail!\n"); exit(10); }
   } else { printf("./a.out <tree_input_file> <taxa_no> <ouput_file>\n"); exit(100); }

   /* --------------- start reading input trees --------------------  */
   k=0; /* keep the number of edges  read */ 
   no_trees=0;

   while (fscanf(tr_input_file, "%s\n", input_nw_str)!=EOF) {
      /* if (no_trees==0) { */
      mix_edge_no=Convert_to_Edges(input_nw_str, taxa_names, I_trees[no_trees], &taxa_no,       I_lf_no);
       Add_Rootedge(mix_edge_no,  I_trees[no_trees], I_lf_no);
       I_trees[no_trees][0][0]=mix_edge_no+1;
      no_trees++;
   }
   fclose(tr_input_file);
   /* --------------- finish reading input trees --------------------  */

   Move_Largest_Tree(no_trees, I_trees); 
   /* move the largest tree to the first so that the first tree is not empty
 *   in the top component */

  for (i=0; i<no_trees; i++)  Print_input_trees(i, I_trees[i][0][0],  I_trees[i]); 
  for (i=0; i<I_lf_no; i++) { printf("%d -->(%s)\n", i+1, taxa_names[i]); }
  
  for (i=0; i<NO_TREES; i++) {
     for (j=0; j<NO_NODES; j++) {
         node_clusters[i][j]=(char *)malloc(I_lf_no+2);
          sorted_tnodes[i][j]=0;
        for (k=0; k<MAX_OUT; k++) I_tree_ch[i][j][k]=0;
        for (k=0; k<NO_LEAVES; k++)  ComponentGraph[j][k]=0;
     }
  }


  /* decomponse subtrees into many disjoint subset of irreducible subtrees */
  /*  trees are indexed from 0,  I_trees[i][0][0] keep the number of edges, 
 *   and edges are in I_trees[i][j], I_lf_no keeps the number of all taxa  */
   Decompose_Trees_Into_Componets(no_trees, I_trees, 
      I_lf_no, node_clusters, sorted_tnodes, I_tree_ch, ComponentGraph);

   no_components=0;
   for (node_ind=2*I_lf_no; node_ind>=1; node_ind--)
       if (ComponentGraph[sorted_tnodes[0][node_ind]][0]>0)  no_components +=1;

   printf("There are %d trees which are decomposed into %d components\n", no_trees,  no_components);

   current_comp=0; final_ntks_list=NULL; 

  /* examine  the tree nodes in the 0-th trees one by one 
 *  the following change is for non-equal size trees           */
  for (node_ind=2*I_lf_no; node_ind>=0; node_ind--){
   /* identify a common component, then work on it */
           /*  printf("here-------Node_ind(%d)---------\n", node_ind); */
   if (ComponentGraph[sorted_tnodes[0][node_ind]][0]>0 || 
   (node_ind==0 &&  Not_All_Empty(no_trees, I_trees)==1)){

   current_comp +=1;
  /* printf("inside_com(%d)  (%d)------\n", node_ind, 
 *  sorted_tnodes[0][node_ind]); */
   
     no_lfs_top=0;
     for (j=0; j<no_trees; j++){
       /* construct the subtrees rooted at the node of the paritition 
         node_clusters[0][sorted_tnodes[0][node_ind]] */
        if (node_ind >0) {
           Construct_Subtrees_Jan2024(I_lf_no, 2*I_lf_no, subtrees[j], 
             I_tree_ch[j], sorted_tnodes[j], node_clusters[j], node_clusters[0],
              sorted_tnodes[0][node_ind], ComponentGraph); 

           /*  Print_input_trees(j, subtrees[j][0][0], subtrees[j]); */ 

     /* Remove all subtrees in common comonpents so that we can deal with the 
      top subtrees which may not have the same number of leaves */
            Remove_Subtree_Edges(I_trees[j],  subtrees[j]);
        /* Print_input_trees(j, I_trees[j][0][0], I_trees[j]); */
        } else {
           Copy_Subtrees_Top_Jan2024(I_trees[j], subtrees[j]);
        }
        /* this part is added for trees of unequal size */
        for (kkk=0; kkk<=I_lf_no; kkk++) { leaf_mapping[j][kkk]=-1; }
        Simplify_Tree_Initialization_Jan2024(subtrees[j],leaf_mapping[j], 
          &no_lfs_top,  top_clusters, node_clusters[j]);
        for (kkk=0; kkk<=I_lf_no; kkk++) { 
          if ( leaf_mapping[0][kkk]<0 && leaf_mapping[j][kkk]>0){ 
                leaf_mapping[0][kkk]=leaf_mapping[j][kkk];
          }
        }
        /*
        for (kkk=0; kkk<=I_lf_no; kkk++) {
          printf("--leaf-mapping--%d  %d\n", kkk,  leaf_mapping[0][kkk]);
        }
        printf("\n");
        */
     } 


     for (j=0; j<no_trees; j++) {
          t_p[j]=(int *)malloc(sizeof(int)*(2*I_lf_no+2));
          for (kkk=0; kkk<2*I_lf_no+2; kkk++) t_p[j][kkk]=-1;
          Process_Tree_Jan2024_v2(subtrees[j], no_lfs_top,  t_ch[j], t_p[j]);
          Print_input_trees(j, subtrees[j][0][0], subtrees[j]);

          subtr_lf_no=no_lfs_top;
         taxa_inci[j]=(int *)malloc(sizeof(int)*(2*subtr_lf_no+1));
         taxa_incidence_str[j]=(char **)malloc(sizeof(char *)*(subtr_lf_no+1));
         for (i=0; i<=subtr_lf_no; i++) {
           taxa_incidence_str[j][i]=(char *)malloc(subtr_lf_no+1);
           taxa_inci_int[j][i]=(int *)malloc(sizeof(int)*(subtr_lf_no+1));
         }

    }

    printf("COMPONENT(%d) lf_no(%d)\n", current_comp, subtr_lf_no);
      
      /* test purpose */
      /*
      if (node_ind==0 ) { printf("here\n"); exit(100); }
      */


      for (j=0; j<subtr_lf_no; j++) zh[j]='0'+j; 
      zh[subtr_lf_no]='\0';

      Q0_front=NULL; Q0_end=NULL;  Q1_front=NULL; Q1_end=NULL;
      kkk=0; it_times=subtr_lf_no/2;

      /* how record the permutations with a total lenght */
      ntks_counts=(int *)malloc(sizeof(int)*(no_trees)*subtr_lf_no);

   
     
       
     for (kkk=0; kkk < it_times; kkk++ ) {
         if (kkk==0) { 
             Qinsert_23(&Q0_front, &Q0_end, zh); 
             cut_off=(no_trees-1)*subtr_lf_no;
         } else { Q0_front=Q_front; Q0_end=Q_end; }


         while (Qdelete(&Q0_front, &Q0_end, &ret_val,  zh)==1) {
          
             if (ret_val <= cut_off) {
                for (k1=2*kkk; k1<subtr_lf_no; k1++) {
                   Generate_Permut(zh, 2*kkk, k1, new_zh);
                   Qinsert_23(&Q1_front, &Q1_end, new_zh);
                }
             } /* if */
         }



         if (kkk !=it_times-1 || subtr_lf_no%2!=0) {
           while (Qdelete(&Q1_front, &Q1_end, &ret_val,  zh)==1) {
             for (k1=(2*kkk+1); k1<subtr_lf_no; k1++) {
                Generate_Permut(zh, 2*kkk+1, k1, new_zh);
                Qinsert_23(&Q0_front, &Q0_end, new_zh);
             }
           }
        } /* if condition */

      

        Q_front=NULL; Q_end=NULL;
        if (kkk !=it_times-1 || subtr_lf_no%2!=0) { 
           Q_front=Q0_front; Q_end=Q0_end; pos=2*kkk+1;
       } else { 
            pos=2*kkk; Q_front=Q1_front; Q_end=Q1_end;
       }

       
       if (Q_front!=Q_end) { 
         printf("      Count in the middle pos(%d)\n", pos); 
         Count_NNN(Q_front, Q_end); 
       }

       q_flag=0; ctt=0; 
       Qptr=Q_front;
     
       for (j=0; j<(no_trees)*subtr_lf_no; j++) ntks_counts[j]=0;   
       

       while (Qptr != NULL) {
          ctt++;

          strcpy(zh, Qptr->item);
          for (j=1; j<=subtr_lf_no; j++){  permut[j]=1+(zh[j-1]-'0'); }
          for (j=0; j<no_trees; j++) {
           Compute_Incidence_July30(permut, taxa_inci[j],  t_p[j], subtr_lf_no);
           Convert_Str_Format_July30(taxa_incidence_str[j], taxa_inci[j],
                subtr_lf_no, t_p[j]);

          } /* j loop for subtrees */


         current_ret_no=Search_minNTKS_Round1(no_trees, subtr_lf_no, permut,
          taxa_incidence_str, pos+1, 0);


         Qptr->score=current_ret_no;
         ntks_counts[current_ret_no]=1+ntks_counts[current_ret_no];
         if (q_flag==0) {
           max=current_ret_no; min_ret_no= current_ret_no; q_flag=1;
         } else {
           if (max< current_ret_no) max=current_ret_no;
           if (min_ret_no > current_ret_no) min_ret_no= current_ret_no;
         }
         Qptr=Qptr->next;
       } /* Qptr loop */

  
       cut_off=Comput_Cutoff(ntks_counts, min_ret_no, max, ctt);
       cut_off=Min(min_ret_no+8, cut_off); 
       if (ctt>CUT_OFF_SIZE) {
          ctt=ctt/(2*CUT_OFF_SIZE); 
          if (ctt==0) ctt=1;
          Remove_BadPermut(&Q_front, &Q_end, cut_off, min_ret_no, ctt);
       }
      
      /* Count_NNN(Q_front, Q_end); */
   } /* kkk loop */

       free(ntks_counts); ntks_counts=NULL;
       

   
    printf("     After sampling search--xxx- min_ret_no(%d)\n", min_ret_no);
    /*  if (Q_front!=Q_end) {  Count_NNN(Q_front, Q_end); } */


   count_ntks=0;
   Qptr=Q_front; new_min_ret_no=min_ret_no;

   q_flag=0;
   while (Qptr!=Q_end) {  /* ---Q */
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
    Qinsert_23(&(Second_Q[0][new_min_ret_no]), &(Second_Q[1][new_min_ret_no]), zh);
        } else if (new_min_ret_no == current_ret_no) {
          
          count_ntks++;
    Qinsert_23(&(Second_Q[0][new_min_ret_no]), &(Second_Q[1][new_min_ret_no]), zh);
        }
    } /* if */
     Qptr=Qptr->next;
   } /* ---  Q while */ 


      min_ret_no=new_min_ret_no;
      printf("       The current networks with ret_no(%d) have been obtained after round 2\n", new_min_ret_no);
     /*  Count_NNN(Second_Q[0][min_ret_no], Second_Q[1][min_ret_no]); */ 



     ntk_list=NULL;
     // gph[0]=(int *)malloc(sizeof(int )*(3*subtr_lf_no+2*min_ret_no));
     // gph[1]=(int *)malloc(sizeof(int )*(3*subtr_lf_no+2*min_ret_no));
      gph[0]=(int *)malloc(sizeof(int )*MAXSIZE);
      gph[1]=(int *)malloc(sizeof(int )*MAXSIZE);

     
      xxx=0;
     while (Qdelete(&(Second_Q[0][min_ret_no]),&(Second_Q[1][min_ret_no]), &ret_val, zh)==1){
        if (xxx>100) break;
      /* printf("------permuation(%d)-- (%s)\n", xxx, zh); */
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

     /*  printf("------xxx_value(%d)--edges(%d)--nodes(%d)\n", xxx, no_gph_edges, no_ntk_nodes); */
      Simplify_Ntks_July22(gph, no_gph_edges);
      /*  Print_Ntk(gph, no_gph_edges); */
      ntk_list=Insert_v2(ntk_list, gph, no_gph_edges, no_ntk_nodes,
           no_gph_edges-no_ntk_nodes+1, subtr_lf_no, zh);
      xxx++;
    } /* while */

  
    free(gph[0]); gph[0]=NULL;
    free(gph[1]); gph[1]=NULL;
    total_rets=total_rets+no_gph_edges-no_ntk_nodes+1;

    final_ntks_list=Combine_Solutions_On_Parts(final_ntks_list, 
             ntk_list, subtr_lf_no, leaf_mapping[0], current_comp, node_ind);


     free(Second_Q[0]); Second_Q[0]=NULL; free(Second_Q[1]); Second_Q[1]=NULL;

    for (j=0; j<no_trees; j++) { 
      for (i=0; i<subtr_lf_no+1; i++){
         if (taxa_incidence_str[j][i]!=NULL) free(taxa_incidence_str[j][i]);  
         if (taxa_inci_int[j][i]!=NULL) free(taxa_inci_int[j][i]); 
     }
     free(t_p[j]);
     free(taxa_inci[j]);
      //free(taxa_inci_int[j]);
    }
   } /* component */
 }    /* node */
  

 /* print out multiple networks, whose size information can be foudn in 
 *  Progress_Log */

 Print_Ntks_500(final_ntks_list, total_rets, taxa_names, I_lf_no, ntk_file);

 for (i=0; i<NO_TREES; i++) {
     for (j=0; j<NO_NODES; j++) {
         free(node_clusters[i][j]);
     }
 }
 fclose(ntk_file);

} /* end main */
