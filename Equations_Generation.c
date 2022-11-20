/* Copyright@Louxin Zhang
 * This program takes a set of weighted trees and a tree-child network that contains 
 * every tree and has  the minimum  number of reticulations, among all tree-child networks
 * containing all the trees and outputs a set of equations for computing the branch 
 * weights to fit the trees in terms of the least squares method. 
 *
 * The file containing input trees list the trees one by one. For each tree, all edges are
 * listed with one edge per row. In each row, the tail and head of each edge and weight 
 * are separated by a space " ", such as "1 12 0.5" represents an edge from 1 to 12 of 
 * weight 0.5. 
 *
 * We also assume that trees nodes are represented by integers. In particular, the root 
 * is 0, which has out degree 1 and the leaves are numbered from 1 to k if there are
 * k leaves.
 * Other nodes are numbered with integeres > k. Since each tree has the same number of 
 * edges, which is 2 times the no. of leaves minus 1, including the edges from 0. Edges of 
 * trees are listed without any separete lines. 
 * For example, a weighted tree with  root 0 , there leaves 1, 2, 3 and internal nodes 7, 8 is 
 * listed as:
 * 0 8 1.00
 * 8 1 0.4
 * 8 7 0.5
 * 7 3 0.9
 * 7 2 0.2
 *
 * Similarly, the edges of a unweigthed binark networks is listed one-by-one.
 *
 * The output file has two parts: 
 *  -- The first is the list of network edges indexed from Edge 1 to the last edges
 *  -- The seond part is a set of equestions which are each derived from an edge weight
 *
 * 
 *
 * Compiling command:  gcc  Equation_generation.c  -o Weight_Eq 
 * Run command:
 *  ./Weight_Eq <tree_input_file> <ntk_input_file> no_leaves no_ret <output_file> 
 */



#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define TREE_OR_LEAF 0
#define RETNODE 1
#define T_EDGE 'a'
#define RED_EDGE 'b'
#define T_LEAF -1
#define T_ROOT 2
#define T_OTHER 1

#define MAX_DEG 3
#define MAX_NO_TREES 300
#define NO_NODES 100 
#define MAXSIZE  100   /* deal with ntks with at most 20 leaves */

typedef struct {
    short lr_ends; /*leftx100+right */
} edge;


typedef struct graph {
    edge *gph;
    struct graph *next;
} g_set;

typedef struct compressed_gph {
    short tree_id;
    short  *extend;  /*  tr*100+ good_tree */
    struct compressed_gph *next;
} zip_gphs;






void  Print_graph3_1(short graph_cpy[][2], short  new_size){
    short i;
    for (i=1; i<=new_size; i++){
      printf("%d %d\n", graph_cpy[i][0], graph_cpy[i][1]);
    }
}

void top_sort(short graph[][2], short node_no, short ntk_edge_no, 
   short sorted_nodes[], short ch[][MAX_DEG], short p[][MAX_DEG], 
   short indicator[]){
    short i, j, k;
    short  indegree[NO_NODES];
    short out_d[NO_NODES];
    short  child[NO_NODES][MAX_DEG]; /* 0 contains the number of children */
    short  roots[NO_NODES];
    short  r_no, no, u;

    for (i=0; i<=node_no; i++) {
            indegree[i]=0;out_d[i]=0;ch[i][0]=0; p[i][0]=0;
    }


    for (i=1; i<=ntk_edge_no; i++) {
        indegree[graph[i][1]] +=1;
        child[graph[i][0]][out_d[graph[i][0]]]=graph[i][1]; 
            /* kth child is child[i][k] */
        out_d[graph[i][0]] +=1; 
        ch[graph[i][0]][0] +=1; p[graph[i][1]][0] +=1;
        k=ch[graph[i][0]][0]; ch[graph[i][0]][k]=graph[i][1]; 
        k=p[graph[i][1]][0]; p[graph[i][1]][k]=graph[i][0]; 
    }

    r_no=0;
    for  (i=0; i<=node_no; i++) { 
       if (indegree[i]==0 &&  out_d[i]==1) 
           {
           roots[r_no]=i; r_no=1; 
           indicator[i]=4; 
           }  /* root */
       else if (indegree[i]==1 &&  out_d[i]==2) indicator[i]=3; /* tree node */
       else  if (indegree[i]==1 &&  out_d[i]==1) indicator[i]=2; /* degree 2*/
       else if (indegree[i]==1 &&  out_d[i]==0) indicator[i]=1;  /* leaf */
       else if (indegree[i]==2 &&  out_d[i]==1) indicator[i]=0;  /* ret */
       else indicator[i]=-1;
    }

    
    no=0;
    while (r_no >0) {
      u=roots[r_no-1]; sorted_nodes[no]=u; 
      no +=1; r_no=r_no-1;
      for (i=0; i<out_d[u]; i++) {
         indegree[child[u][i]]=indegree[child[u][i]]-1;
         if (indegree[child[u][i]]==0) {roots[r_no]=child[u][i]; r_no +=1; }
      }
    } /* whilde loop */
}

void Procese100(short tree[][2], short edges_no, short nodes_no, short ch[][MAX_DEG], 
   short p[][MAX_DEG]){
    short  i, j, k;

   for (i=0; i<=nodes_no;  i++) { ch[i][0]=0; p[i][0]=0; }
   for (i=0; i<edges_no;  i++) { 
      ch[tree[i][0]][0] +=1; k=ch[tree[i][0]][0]; ch[tree[i][0]][k]=tree[i][1];
      p[tree[i][1]][0] +=1; k=p[tree[i][1]][0]; p[tree[i][1]][k]=tree[i][0];
   }
} /* process100 */

short max(short x, short y){ if (x>y) return x; else return y; }

short  Check_On_Tree(short lf1,  short lf2, short ch[][MAX_DEG], short t_p[][MAX_DEG],
     short  i, short node_map[]) {
     short p,p1;

     p=t_p[lf2][1]; 
     if (ch[p][1]==lf1 || ch[p][2]==lf1) {
       node_map[p]=i;
       p1=t_p[p][1]; 
       if (ch[p1][1]==p) ch[p1][1]=max(lf1, lf2); else ch[p1][2]=max(lf1, lf2);
       if (lf1 > lf2) t_p[lf1][1]=p1; else t_p[lf2][1]=p1;
       return 1;
     } else return 0;
} /* check on tree */

void Printf_Tree_Image(FILE *out, short image[][2], short image_edge_no){
   short i;
   fprintf(out, "Tree image\n");
   for (i=0; i<image_edge_no; i++) fprintf(out, "%d %d\n", image[i][0], image[i][1]);
   fprintf(out, "\n\n");
}

/*
short Print_node_image(FILE *out, short node_image[], short a){
    short i;

   fprintf(out, "Tree node mapping\n");
   for (i=0; i<a; i++) 
     if (node_image[i]!=-1) fprintf(out, "%d --> %d\n", i, node_image[i]);
   fprintf(out, "\n\n");
}
*/

/* edge index start from 1 */
short Comput_Edge_Index(short ntk[][2], short no, short p, short son){
    short i; 
    
    for (i=1; i<=no; i++)  { if (p==ntk[i][0] && son==ntk[i][1]) return i; }
}


short  Tree_Containment(short t_edges_no, short t_child[][MAX_DEG], 
     short t_parent[][MAX_DEG], short ntk_ch[][MAX_DEG], short ntk_p[][MAX_DEG], 
     short sorted[], short node_no, short node_indicator[], short image[][2],
     short *i_edge_no, short node_map[]){
     short i, j, k;
     short c1, c2, p2, p3;
     short ch[NO_NODES][MAX_DEG], p[NO_NODES][MAX_DEG];/* ntk copy */
     short t_ch[NO_NODES][MAX_DEG], t_p[NO_NODES][MAX_DEG];
     int t_no, ind[NO_NODES], image_edge_no;

     
     for (i=0; i<=node_no; i++) {
         ch[i][0]=ntk_ch[i][0];
         for (j=1; j<=ch[i][0]; j++) { ch[i][j]=ntk_ch[i][j]; }
         p[i][0]=ntk_p[i][0];
         for (j=1; j<=p[i][0]; j++) { p[i][j]=ntk_p[i][j]; }

         ind[i]=node_indicator[i];
  
         /* tree_part */
         t_ch[i][0]=t_child[i][0];
         for (j=1; j<=t_ch[i][0]; j++) { t_ch[i][j]=t_child[i][j]; }
         t_p[i][0]=t_parent[i][0];
         for (j=1; j<=t_p[i][0]; j++) { t_p[i][j]=t_parent[i][j]; }
     }

     t_no=t_edges_no/2;  /* the number of internal nodes */
     image_edge_no=0;
     for (i=0; i<=node_no; i++) node_map[i]=-1;
   
      /* printf("t_no %d\n", t_no);  */
     for (j=node_no-1; j>=0; j--) {
        i=sorted[j]; 
        switch (ind[i]) {
          case 4:
             if (t_no==1) return 1; else return 0;
             break;
          case 3:
             /* printf("%d: tree node %d %d\n", i, ch[i][1], ch[i][2]); */
             p3=p[i][1]; if (ch[p3][1]==i) k=1; else k=2;
             if (ind[ch[i][1]]==0 || ind[ch[i][2]]==0) {
                   /* printf("   1 ret 1 leaf\n"); */  
                if (ind[ch[i][1]]==0) { c1=ch[i][1]; c2=ch[i][2];
                } else {c2=ch[i][1]; c1=ch[i][2];}
                if (p[c1][1]!=i) p2=p[c1][1]; else p2=p[c1][2];
                     /* printf("p2(%d) p3(%d)\n", p2, p3); */
                if (Check_On_Tree(c2, ch[c1][1], t_ch, t_p, i, node_map)==1) {
                        /*  printf("           merge\n");  */ 
                      /* image update */ 
                   image[image_edge_no][0]=p3; image[image_edge_no][1]=i;
                   image_edge_no +=1;
                   image[image_edge_no][0]=i; image[image_edge_no][1]=c1;
                   image_edge_no +=1;
                         
                   ch[p3][k]=max(ch[c1][1],c2); 
                     /* modify the degre of p2 and child list */
                   if (ch[p2][1]==c1) { ch[p2][1]=ch[p2][2];}
                   ch[p2][0]=1; ind[p2]=2; /* p2 becom degree 2 */
                     /* modify the i into a leaf */
                   if (ch[c1][1] >= c2) p[ch[c1][1]][1]=p3; else p[c2][1]=p3;
                    t_no=t_no-1;
                       /* printf(" t_no: %d\n", t_no);  */
                   if (t_no==0) {
                      *i_edge_no=image_edge_no;
                      /*
                       Printf_Tree_Image(image, image_edge_no);
                       Print_node_imge(node_map, NO_NODES);
                       */
                       return 1; 
                      }
                  }  else { /* child of c1 goes to c2 */
                      /* printf("           split\n"); */  

                      image[image_edge_no][0]=p3; image[image_edge_no][1]=i;
                      image_edge_no +=1;
                      image[image_edge_no][0]=p2; image[image_edge_no][1]=c1;
                      image_edge_no +=1;
                      
                     if (ch[p2][1]==c1) { ch[p2][1]=ch[c1][1]; 
                     } else { ch[p2][2]=ch[c1][1]; }
                      p[ch[c1][1]][1]=p2;
                      
                     ch[p3][k]=c2; p[c2][1]=p3;
                  }
             } else if (ind[ch[i][1]]==1 &&  ind[ch[i][2]]==1) {
                      /* printf("     two leaves\n");  */ 
                  /* both are leaves */
                  c1=ch[i][1]; c2=ch[i][2];
                  
                  if (Check_On_Tree(c2, c1, t_ch, t_p, i, node_map)==1) {
                      image[image_edge_no][0]=p3; image[image_edge_no][1]=i;
                      image_edge_no +=1;

                      ch[p3][k]=max(c1, c2);
                      if (c1>c2) p[c1][1]=p3; else p[c2][1]=p3;

                      t_no=t_no-1; 
                      if (t_no==0) {
                          *i_edge_no=image_edge_no;
                          /*
                          Printf_Tree_Image(image, image_edge_no);
                          Print_node_imge(node_map, NO_NODES);
                           */
                          return 1; 
                      }
                  } else return 0;
                  
             } 
             break;
          case 2: /* degree 2 */
             p3=p[i][1];
             image[image_edge_no][0]=p[i][1]; image[image_edge_no][1]=i;
             image_edge_no +=1;
              
             if (ch[p3][1]==i) ch[p3][1]=ch[i][1]; else ch[p3][2]=ch[i][1];
             break;
          case 1: /* leaf */  
             image[image_edge_no][0]=p[i][1]; image[image_edge_no][1]=i;
             image_edge_no +=1; node_map[i]=i;
             break;
          case 0: /* ret */
              break;
          default: 
             printf("unusual case! %d\n", i); 
             exit(100);
        }
     }
}



void  Reverse_CP(edge *g, short graph[][2], short edge_no){
     short i;

     for (i=0; i<edge_no; i++) {
         graph[i+1][0]=g[i].lr_ends/100; graph[i+1][1]=g[i].lr_ends%100;
     }
}



void   Process3_1(short graph[][2], short ret_edges[][2], short tree_edges[][2],
  short leaves[],  short edges, short node_no){
  short i, j, r_ind, t_ind, lf_ind;
  short node_indicator[NO_NODES], indicator1[NO_NODES];
  short temp[2];

  for (i=0; i<=node_no; i++) { node_indicator[i]=-1; indicator1[i]=0; }
        node_indicator[0]=0;
  for (i=1; i<=edges; i++) {
    node_indicator[graph[i][1]] =node_indicator[graph[i][1]]+1;
    indicator1[graph[i][0]] =indicator1[graph[i][0]]+1;
  }

  lf_ind=0;
  for (i=0; i<=node_no; i++) {
     if (indicator1[i]==0 && node_indicator[i]==TREE_OR_LEAF)
         {  leaves[lf_ind]=i; lf_ind +=1; }
  }


  r_ind=0; t_ind=0;
  for (i=1; i<=edges; i++) {
     if  (node_indicator[graph[i][1]]==TREE_OR_LEAF){
        tree_edges[t_ind][0]=graph[i][0]; tree_edges[t_ind][1]=graph[i][1];
        t_ind=1+ t_ind;
     }
     else  if (node_indicator[graph[i][1]]==RETNODE){
        ret_edges[r_ind][0]=graph[i][0]; ret_edges[r_ind][1]=graph[i][1];
        r_ind=1+ r_ind;
     }
  }

  for  (i=0; i<r_ind; i++) {
      for (j=i+1; j<r_ind; j++) {
          if (ret_edges[j][1]==ret_edges[i][1]){
             temp[0]=ret_edges[j][0]; temp[1]=ret_edges[j][1];
             ret_edges[j][0]=ret_edges[i+1][0]; ret_edges[j][1]=ret_edges[i+1][1];
             ret_edges[i+1][0]=temp[0]; ret_edges[i+1][1]=temp[1];
             break;
          }
      }
  }
 } /* process3 */


short  check_edge(short a, short b,  short edges[][2], short no_edges){
      short i, j;

  for (i=0; i<no_edges; i++) { if (a==edges[i][0] && b==edges[i][1]) { return 1;} }
      return 0;
}

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
    flag = 1+flag;
    if (flag > 50) exit(100);
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


} /* end */



short  Check_Iso(short gph1[][2], g_set *ntks, short no_eg, short node_no, short ret_no,  short lf_no){
     short  map[NO_NODES];
     g_set *ntk_ptr;
      short gph2[MAXSIZE][2];
      short ret_egs1[MAXSIZE][2], ret_egs2[MAXSIZE][2];
      short tree_egs1[MAXSIZE][2], tree_egs2[MAXSIZE][2];
      short  array_lf1[NO_NODES], array_lf2[NO_NODES];
      short i;


     for (i=0; i<=node_no; i++) map[i]=-1;

     Process3_1(gph1, ret_egs1, tree_egs1, array_lf1, no_eg, node_no);

     ntk_ptr=ntks;
     while (ntk_ptr!=NULL) {
       Reverse_CP(ntk_ptr->gph, gph2, no_eg);
       /* Print_graph3_1(gph2, no_eg); */
       Process3_1(gph2, ret_egs2, tree_egs2, array_lf2, no_eg, node_no);
       /*  Make_Map(map, lf_no);  */
       if (Is_Isomorphism(map, ret_egs1, tree_egs1, ret_egs2, tree_egs2,
         array_lf1, node_no, ret_no, lf_no)==1) { 
            /* ntk_ptr->freq +=1; */ return 1; 
        }
       ntk_ptr=ntk_ptr->next;
     }
     return 0;
} 




void Print_Coe_Matrix_WT(FILE *out_file, short tree_edges_no, short  ntk_edge_no,
   short coe_matrix[][MAXSIZE], float wt[], short ntk_p[][MAX_DEG], short graph[][2]){
   short i, j;

     /* fprintf(out_file, "Coe Matrix\n"); */

   for (i=1; i<=tree_edges_no; i++){
     for (j=1; j<=ntk_edge_no;  j++) { 
       if (ntk_p[graph[j][0]][0]!=2)
        fprintf(out_file, "%d ", coe_matrix[i][j]); 
     }
     fprintf(out_file, "  %1.2f\n", wt[i]);
   }
}

void  Define_Equations(short ntk[][2], short ntk_edge_no, short image[][2],
      short image_edge_no, short tree[][2], short tree_edge_no,
      short node_map[], short coefficient[][MAXSIZE], short ntk_p[][MAX_DEG]){
      /* tree_edge_no by ntk_edge_no */

      short i, j;
      short image_p[NO_NODES];   /* parent relation in the image tree */
      short  a, b, x, y;
      short index;

      for (i=0; i<NO_NODES; i++) { image_p[i]=-1; }
      for (i=0; i<image_edge_no; i++) { image_p[image[i][1]]=image[i][0]; }
      for (i=0; i<=tree_edge_no; i++) { 
            for (j=0; j<=ntk_edge_no; j++) coefficient[i][j]=0;
      }


      for (i=1; i<=tree_edge_no; i++) {
          a=node_map[tree[i][0]]; b=node_map[tree[i][1]]; x=b;
          while (x != a) {
             y=image_p[x]; 
             index=Comput_Edge_Index(ntk, ntk_edge_no, y, x);
             if (ntk_p[y][0]!=2) coefficient[i][index]=1; x=y;
          }
      }
} /* define equations */


void Display_Trees_CoeMatrix_WT(short trees[][NO_NODES][2], float trees_wt[][NO_NODES],
  short no_trees, short t_edges_no, short t_ch[][NO_NODES][MAX_DEG],
    short t_p[][NO_NODES][MAX_DEG],
     short graph[][2], short ntk_edge_no, short ret_no,  FILE *outfile){

   short i,j;
   /*
   short graph[MAXSIZE][2];
    */
   short node_no;
   short ntk_ch[NO_NODES][MAX_DEG], ntk_p[NO_NODES][MAX_DEG];
   short indicator[NO_NODES], node_ind[NO_NODES];
   short  sorted[NO_NODES];
   short image[MAXSIZE][2];
   short i_edge_no;
   short tr_ntk_node_map[NO_NODES];
   short  coe_matrix[NO_NODES][MAXSIZE];




      /*
      Reverse_CP(gph, graph, ntk_edge_no);
      Print_graph3_1(graph, ntk_edge_no); printf("inside =====111=====\n");
      */
       node_no=ntk_edge_no-ret_no+1;   /* include node 0*/

       top_sort(graph, node_no, ntk_edge_no, sorted, ntk_ch, ntk_p, indicator);

        for (j=0; j<node_no; j++) {
            printf("%d ", sorted[j]); 
        }
        printf("\n");

            fprintf(outfile, "network:\n");
        for (j=1; j<=ntk_edge_no; j++)
            fprintf(outfile,"edge(%2d): %2d  %2d\n", j,  graph[j][0],graph[j][1]);

        for (i=0; i<no_trees; i++) {
          if (Tree_Containment(t_edges_no, t_ch[i], t_p[i], ntk_ch, ntk_p,
        sorted, node_no, indicator, image, &i_edge_no, tr_ntk_node_map )==1){
            printf("-------------Process Tree %d\n", i);    
            Define_Equations(graph, ntk_edge_no, image, i_edge_no,
             trees[i], t_edges_no, tr_ntk_node_map, coe_matrix, ntk_p);
            fprintf(outfile, "Equations from Tree %d\n", i);
         Print_Coe_Matrix_WT(outfile, t_edges_no, ntk_edge_no, coe_matrix, trees_wt[i], ntk_p,
            graph);
            fprintf(outfile, "\n");
        }
       } /* for loop on trees */
            fprintf(outfile, "\n");
    
} /* weightvertion */

short  Simplify_Ntks(short gph[][2], short no_edges, short node_no){
    short i, j, k, r;
    short  indegree[NO_NODES],out_d[NO_NODES],internal[NO_NODES],external[NO_NODES];
    short  v_type[NO_NODES], cpy[2];
    short no_lfs, no_trs_rets;
    float wt_cp;

    for (i=0; i<NO_NODES; i++) { indegree[i]=0; out_d[i]=0; }
    for (i=1; i<=no_edges; i++) { indegree[gph[i][1]] +=1; out_d[gph[i][0]] +=1;}
    k=0; j=1;
    for  (i=0; i<NO_NODES; i++) {
       if (indegree[i]==0 &&  out_d[i]==1) { v_type[i]=4; }  /* root */
       else if (indegree[i]==1 &&  out_d[i]==2) {
            v_type[i]=3; internal[i]=k; k +=1;
       } /* tree node */
       else  if (indegree[i]==1 &&  out_d[i]==1) {
                v_type[i]=2; /* degree 2*/
                internal[i]=k; k +=1;
       } else if (indegree[i]==1 &&  out_d[i]==0) {
               v_type[i]=1;  /* leaf */
               external[i]=j; j +=1;
       } else if (indegree[i]==2 &&  out_d[i]==1) {
                v_type[i]=0;  /* ret */
                internal[i]=k; k +=1;
       } else v_type[i]=-1;
    }
   no_lfs=j+1; no_trs_rets=k;
   for (i=1; i<=no_edges; i++) {
          if (v_type[gph[i][0]]==3 ||v_type[gph[i][0]]==0)
              gph[i][0]= no_lfs + internal[gph[i][0]];
          if (v_type[gph[i][1]]==3 || v_type[gph[i][1]]==0)
              gph[i][1]= no_lfs + internal[gph[i][1]];
          if (v_type[gph[i][0]]==4) { r=i; gph[i][0]=0; }
          if (v_type[gph[i][1]]==1) gph[i][1]=external[gph[i][1]];
    }
    cpy[0]=gph[1][0]; cpy[1]=gph[1][1];
    gph[1][0]=gph[r][0]; gph[1][1]=gph[r][1];
    gph[r][0]=cpy[0]; gph[r][1]=cpy[1];
 
    gph[0][0]=gph[no_edges][0]; gph[0][1]=gph[no_edges][1];


} /* unweighted version */



short  Simplify_Ntks_WT(short gph[][2], short no_edges, short node_no, 
     float wt[]){
    short i, j, k, r;
    short  indegree[NO_NODES],out_d[NO_NODES],internal[NO_NODES],external[NO_NODES];
    short  v_type[NO_NODES], cpy[2];
    short no_lfs, no_trs_rets;
    float wt_cp;

    for (i=0; i<NO_NODES; i++) { indegree[i]=0; out_d[i]=0; }
    for (i=1; i<=no_edges; i++) { indegree[gph[i][1]] +=1; out_d[gph[i][0]] +=1;}
    k=0; j=1;
    for  (i=0; i<NO_NODES; i++) {
       if (indegree[i]==0 &&  out_d[i]==1) { v_type[i]=4; }  /* root */
       else if (indegree[i]==1 &&  out_d[i]==2) {
            v_type[i]=3; internal[i]=k; k +=1;
       } /* tree node */
       else  if (indegree[i]==1 &&  out_d[i]==1) {
                v_type[i]=2; /* degree 2*/
                internal[i]=k; k +=1;
       } else if (indegree[i]==1 &&  out_d[i]==0) {
               v_type[i]=1;  /* leaf */
               external[i]=j; j +=1;
       } else if (indegree[i]==2 &&  out_d[i]==1) {
                v_type[i]=0;  /* ret */
                internal[i]=k; k +=1;
       } else v_type[i]=-1;
    }
   no_lfs=j+1; no_trs_rets=k;
   for (i=1; i<=no_edges; i++) {
          if (v_type[gph[i][0]]==3 ||v_type[gph[i][0]]==0)
              gph[i][0]= no_lfs + internal[gph[i][0]];
          if (v_type[gph[i][1]]==3 || v_type[gph[i][1]]==0)
              gph[i][1]= no_lfs + internal[gph[i][1]];
          if (v_type[gph[i][0]]==4) { r=i; gph[i][0]=0; }
          if (v_type[gph[i][1]]==1) gph[i][1]=external[gph[i][1]];
    }
    cpy[0]=gph[1][0]; cpy[1]=gph[1][1];
    wt_cp=wt[1];
    gph[1][0]=gph[r][0]; gph[1][1]=gph[r][1];
     wt[1]=wt[r];
    gph[r][0]=cpy[0]; gph[r][1]=cpy[1];
     wt[r]=wt_cp;
 
    gph[0][0]=gph[no_edges][0]; gph[0][1]=gph[no_edges][1];
    wt[0]=wt[no_edges];


} 
/* simplify_process2 */



void Print_Ntks_Weight(short trees[][NO_NODES][2], float trees_wt[][NO_NODES], 
   short no_trees, short t_ch[][NO_NODES][MAX_DEG],
    short t_p[][NO_NODES][MAX_DEG],
   g_set *new_ntks_cp, short ret_no, short lf_no,  FILE *outfile){
       short graph[MAXSIZE][2];
       short no_eg, node_no;
       short j;
       g_set *ptr;
       int ttt;

       ptr=new_ntks_cp; 
       no_eg=3*ret_no+2*lf_no-1;
       node_no=2*ret_no+2*lf_no;
    
       ttt=0;

       while (ptr!=NULL) {
          Reverse_CP(ptr->gph, graph, no_eg);
          Print_graph3_1(graph, no_eg); 
          printf("=====before simply =====\n");
          Simplify_Ntks(graph, no_eg,  node_no);
          Print_graph3_1(graph, no_eg); 
          printf("=====after  simply =====\n");
          if (Check_Iso(graph, ptr->next, no_eg, node_no, ret_no, lf_no)==0){
            /* printf("------hhhh-----zz----\n"); */
            ttt +=1;
     Display_Trees_CoeMatrix_WT(trees, trees_wt,  no_trees, 2*lf_no-1,  t_ch,
               t_p, graph, no_eg, ret_no, outfile);
            /* printf("------hhhhxhhhhh-zz--------\n"); */
          }
          ptr=ptr->next;
       }
       printf(" Output %d olution networks\n", ttt);
} /*---*/

void Print_TREE_sept2021(short tree[][2], int edge_no, float wt[]){
   int i;

   printf("tree \n"); 
   for (i=0; i<edge_no; i++) 
       { printf("%d %d %2.2f\n", tree[i][0],  tree[i][1], wt[i]); } 
   printf("\n\n");

}

/* this program takes a set of weighted trees and a tc networks models such that the trees are contained in the
* networks,  number leaves in both trees and networks, num of reticulations, outputfile.
* and output a system of linear equations derived from the embedding of tress in the netwokr.
* Hypothesis:  the network is a network with the minimum number of networks among all the networks displaying
* networks
     ./a.out tree_set_input_file, <network_file> <lf_no> <ret_no>  <output_file>

Remark:  the trees must be displayed in the given network

  Examples:
    input tree file that contains two weighted trees:
      0  8 1.0
10  4 1.7
11  2 1.8
9  12 1.6
8  10 1.0
8  5 3.0
9  11 1.0
10  9 1.0
11  1 0.7
12  3 0.9
12  6 0.8
0  8 1.0
10  4 1.7
11  2 1.8
9  12 1.6
8  10 1.0
8  9 1.6
9  5 1.4
10  11 2.0
11  1 0.7
12  3 0.9
12  6 0.8

 and network file:

0 8
8 12
8 13
9 15
10 16
9 11
10 11
11 18
12 9
13 10
12 14
13 14
14 4
15 1
16 5
15 17
16 17
17 2
18 3
18 6
 *               
 *              
 */


void main(int argc, char *argv[]){

short i, j, k;
int node1, node2;
FILE *network_input_file, *tree_input_file, *out_file; 
short ret_no,  ntk_edge_no, node_no, nonleaf_nodes;
short no_trees, lf_no; 
short  trees[MAX_NO_TREES][NO_NODES][2];
float  trees_wt[MAX_NO_TREES][NO_NODES];
short t_ch[MAX_NO_TREES][NO_NODES][MAX_DEG], t_p[MAX_NO_TREES][NO_NODES][MAX_DEG];

zip_gphs *ntks;  
g_set *new_ntks_cp, *my_ntk;
short my_edge;
float wt;



   if (argc>=6) {
     tree_input_file=fopen(argv[1], "r"); 
     lf_no=atoi(argv[3]); ret_no=atoi(argv[4]);
   } else  printf("input format: 5 paraments at least.\n");


   k=0; no_trees=0;
  while (fscanf(tree_input_file, "%d %d %f\n", &node1, &node2, &wt)!=EOF) {
      k=k+1;
      trees[no_trees][k][0]=(short)node1; trees[no_trees][k][1]=(short)node2;
      trees_wt[no_trees][k]=wt;   /* change for real data */
      if (k==(2*lf_no-1)) { k=0;   no_trees +=1;  }
   }
   fclose(tree_input_file);

   printf("tree no: %d, no. edges: %d\n", no_trees, 2*lf_no-1);

   ntks=NULL;
   ntk_edge_no=2*lf_no -1;
   node_no=2*lf_no;
   for (i=0; i<no_trees; i++) {
     Simplify_Ntks_WT(trees[i], ntk_edge_no, node_no,  trees_wt[i]);
     Print_TREE_sept2021(trees[i], ntk_edge_no, trees_wt[i]);
     /*
      */
     Procese100(trees[i], ntk_edge_no, node_no, t_ch[i], t_p[i]);
   }

   ntk_edge_no=3*ret_no+2*lf_no-1; /* add 1 edge entering the root */
   node_no=2*ret_no+2*lf_no; /*2xlf_no nodes including 0 & increases 2 per ret*/
   nonleaf_nodes=2*ret_no+lf_no-1; 
   /* lf_no-1 orginial internal nodes and increase 2 per ret, not including root */
   if (argc>=6) {
     network_input_file=fopen(argv[2], "r");
   } else { printf("input format: 4 paraments at least.\n"); exit(10); }

   new_ntks_cp=NULL;
   my_ntk=(g_set *)malloc(sizeof(g_set));
   my_ntk->gph=(edge *)malloc(ntk_edge_no*sizeof(edge));
   my_ntk->next=new_ntks_cp;


   k=0;
   while (fscanf(network_input_file, "%d %d\n", &node1, &node2)!=EOF) {
      my_edge=(short)node1*100+(short)node2;
      (my_ntk->gph)[k].lr_ends=my_edge;
       k=k+1;
      if (k== ntk_edge_no) { 
        k=0;   
        new_ntks_cp=my_ntk;
        my_ntk=(g_set *)malloc(sizeof(g_set));
        my_ntk->gph=(edge *)malloc(ntk_edge_no*sizeof(edge));
        my_ntk->next=new_ntks_cp;
      }
   }
   fclose(network_input_file);
   

   out_file=fopen(argv[5], "w");

   if (out_file==NULL) {  printf("opening output file fails\n"); exit(100); }
   Print_Ntks_Weight(trees, trees_wt, no_trees, t_ch, t_p, new_ntks_cp, 
              ret_no, lf_no, out_file);

   close(out_file);
} /* end main */
