package COMPONENT_GRAPHER;
/*
 *  COMPONENT-GRAPHER v1.0
 *  
 *  Copyright (C) 2015-2016  Etienne Lord
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.HashMap;
import java.io.*;
import java.util.ArrayList;
import java.util.Stack;

/**
 * Simple node in the graph
 * @author Etienne Lord, Jananan Pathmanathan
 * @since October/November 2015
 */
public class graph {
    
   HashMap<Integer,HashMap<Integer,Boolean>> adjlist=new HashMap<Integer,HashMap<Integer,Boolean>>();
   HashMap<String,Integer> node_to_id=new HashMap<String,Integer>();
   HashMap<Integer,String> id_to_node=new  HashMap<Integer,String>(); 
   HashMap<Integer,Integer> old_id_to_id=new  HashMap<Integer,Integer>(); 
   HashMap<Integer,Integer> id_to_old_id=new  HashMap<Integer,Integer>(); 
   public int total_nodes=0;
   public int total_edges=0;
   public boolean directed=false;
   public static int infinity=9997999;
   /////////////////////////////////////////////////////////////////////////////
   /// Epath
   int[][] epath=new int[1][1]; // epath matrix
   HashMap<Integer,ArrayList<Integer>> epath_extra=new HashMap<Integer,ArrayList<Integer>>(); 
   int epath_next=total_nodes+1;
   ArrayList<ArrayList<Integer>> paths=new ArrayList<ArrayList<Integer>>();
   /////////////////////////////////////////////////////////////////////////////
   /// Class result
   public class results {
       public float total_len3=0;
       public float total_len4=0;
       public float total_loop3=0;
       public float total_loop4=0;
   }
   /////////////////////////////////////////////////////////////////////////////
   // Constructor
   public graph() {};
   
   public String graph_info() {
       String str="";
       str+="vertex: "+this.total_nodes+"\nedges:"+total_edges+"\n";
       return str;
   }
      
   public boolean load_graph(String filename, boolean directed) {
       this.directed=directed;
       this.adjlist.clear();
       this.node_to_id.clear();
       this.id_to_node.clear();
       total_nodes=0;
       try {
           BufferedReader br=new BufferedReader(new FileReader(new File(filename)));
           String data="";
           while(br.ready()) {
               data=br.readLine();
               String[] stri=data.split("\t");
               if (stri.length>1) {
                    int source=getNode(stri[0]);
                    int dest=getNode(stri[1]);
                    addEdge(source,dest);
                    total_edges++;
                   if (!directed)  addEdge(dest,source);
               }
           }
           br.close();
       } catch(Exception e) {
           return false;
       }
       total_nodes=node_to_id.size();
       return true;
   }
   
   //--Need this for the node renumbering
   public int addNode(int node_id) {
       Integer id=old_id_to_id.get(node_id);
       if (id!=null) return id;
       id=this.old_id_to_id.size();
       this.old_id_to_id.put(node_id,id);
       this.id_to_old_id.put(id,node_id);
       return id;
   }
   
   
   public int getNode(String stri) {
       stri=stri.trim();
       Integer node=this.node_to_id.get(stri);
       if (node!=null) return node;
       node=this.node_to_id.size();
       node_to_id.put(stri,node);
       this.id_to_node.put(node,stri);
       return node;
   }
   
   public void addEdge(int source, int dest) {
       HashMap<Integer,Boolean> tmp=new HashMap<Integer,Boolean>();
       if (this.adjlist.containsKey(source)) {
           tmp=this.adjlist.get(source);
       }
       tmp.put(dest, Boolean.TRUE);
       this.adjlist.put(source,tmp);
   }
   
   public ArrayList<results> findAllLoops(boolean verbose) {
       ArrayList<results> result=new  ArrayList<results>();
       for (int i=0; i<this.total_nodes;i++) {
           result.add(findLoops(i, verbose));
           //System.out.println(id_to_node.get(i)+"\t"+r.total_len3+"\t"+r.total_loop3+"\t"+r.total_len4+"\t"+r.total_loop4);
       }
       return(result);
   }
   
    public float[] findTriplets(boolean verbose) {
       float[] result=new float[total_nodes];
        for (int i=0; i<this.total_nodes;i++) {
           result[i]=find_triplet(i,verbose);
           //System.out.println(id_to_node.get(i)+"\t"+result[i]);
       }
        return result;
   }
   
    
    
   public results findLoops(int s, boolean verbose) {
       
       results result=new results();
       for (int v:get_adj(s).keySet()) {
           for (int w:get_adj(v).keySet()) {
               // Loop of and path len 3
                 if (is_valid3(s,v,w)) {
                    result.total_len3++;
                    if (verbose) System.out.print(id_to_node.get(s)+" "+id_to_node.get(v)+" "+id_to_node.get(w));
                     if (get_adj(s).containsKey(w))  { 
                         if (verbose) System.out.println("*"); 
                         result.total_loop3++;
                     } else { 
                         if (verbose) System.out.println("");
                     }
                 }
                 
                 for (int z:get_adj(w).keySet()) {
                   // Loop and paths of len 4
                    if (is_valid4(s,v,w,z)) {
                       if (verbose)  System.out.print(id_to_node.get(s)+" "+id_to_node.get(v)+" "+id_to_node.get(w)+" "+id_to_node.get(z));
                        result.total_len4++; 
                       if (get_adj(s).containsKey(z)) { 
                             if (verbose) System.out.println("*"); 
                            result.total_loop4++;
                        } else { 
                            if (verbose) System.out.println("");
                        }
                    }
                    
               }
           }
       }
       return result;
   }
   
   HashMap<Integer,Boolean> get_adj(int node) {
       HashMap<Integer,Boolean> tmp=this.adjlist.get(node);
       if (tmp==null) return new  HashMap<Integer,Boolean>();
       return tmp;
   }
   
   boolean is_valid3(int s, int v, int w) {
       return (s!=v&&s!=w&&v!=w);
   }
   
    boolean is_valid4(int s, int v, int w, int z) {
       return (s!=v&&s!=w&&s!=z&&v!=w&&v!=z&&w!=z);
   }
    
    boolean is_valid5(int s, int v, int w, int z) {
       return (s!=v&&s!=w&&s!=z&&v!=w&&v!=z&&w!=z);
   }
    
   public Float find_triplet(int s, boolean verbose) {
       
       float total=0;
       ArrayList<Integer> neighbor=new ArrayList<Integer>();
       for (int i:get_adj(s).keySet()) neighbor.add(i);
       int len=neighbor.size();      
       for (int i=0; i<len;i++)
           for (int j=0;j<len;j++) {
               if (i>j&&!get_adj(neighbor.get(i)).containsKey(neighbor.get(j))) {
                   if (verbose) System.out.println(id_to_node.get(neighbor.get(i))+"\t"+id_to_node.get(s)+"\t"+id_to_node.get(neighbor.get(j)));
                   total++;
               }
           }
       //System.out.println(total);
       return total;
   } 
   
  //http://stackoverflow.com/questions/10226251/how-to-find-the-number-of-different-shortest-paths-between-two-vertices-in-dire
   public int[] BFS_count(int s) {
      int[] tmp=new int[total_nodes+1];
      return tmp;
  } 
   
   public int[] DFS(int s,int ignore) {
       int[] tmp=new int[total_nodes+1];
       boolean[] visited=new boolean[total_nodes];
       int total_to_visit=total_nodes;
       if (ignore!=-1) total_to_visit--;
       for (int i=0; i<total_nodes;i++) {
           visited[i]=false;
           tmp[i]=0;
       }
       int index=0;
       Stack<Integer> stack=new Stack<Integer>();
      
       stack.push(s);
       while(total_to_visit>0) {
        // Find next vertice
        if (stack.isEmpty()) {
            for (int i=0; i<total_nodes;i++) {
                if (!visited[i]&&i!=ignore) { stack.push(i);break;} 
            }
        }   
        while(!stack.isEmpty()) {
            int v=stack.pop();
            if (!visited[v]) {
                visited[v]=true;
                total_to_visit--;
                tmp[v]=index;
                for (int i:get_adj(v).keySet()) 
                    if (!visited[i]&&i!=ignore) stack.push(i);
            }
        }
        index++;
       }
       tmp[total_nodes]=index;
       return tmp;
   }
   
   // Special DFS to find if s is and articulation point
   public boolean is_global_articulation_point(int s) {
       // Take any node != s
       int start=(s+2<total_nodes?s+1:s-1);
       int[] tmp=DFS(start,s);
       //System.out.println(id_to_node.get(s)+" "+tmp[total_nodes]);
      if (tmp[total_nodes]>1) return true;
       return false;
   }
   
    // Special DFS to find if s is and articulation point
   public boolean is_local_articulation_point(int s) {
       // Construct local graph for s
       subgraph sub=new subgraph(this,s);
       // get the first node_id..
         int[] tmp=sub.DFS(0,-1);
       return tmp[sub.total_nodes]>1;
     
   }
   
   public ArrayList<ArrayList<Integer>> getCC() {
       ArrayList<ArrayList<Integer>> tmp =new  ArrayList<ArrayList<Integer>>();
      
       int[] data=DFS(0,-1);
        for (int i=0;i<data[total_nodes];i++) tmp.add(new ArrayList<Integer>());        
       for (int i=0; i<data.length-1;i++) {
           int group=data[i];
            ArrayList<Integer>tmp2=tmp.get(group);
           tmp2.add(i);
           tmp.set(group, tmp2);
       }
       
       return tmp;
   }
   
   /**
    * Brandes, U. (2001). A faster algorithm for betweenness centrality*. Journal of Mathematical Sociology, 25(2), 163-177.
    * See also:
    * Brandes, U. (2008). On variants of shortest-path betweenness centrality and their generic computation. Social Networks, 30(2), 136-145.
    * @return array of the betweenness of each node
    */
   public float[] Betweenness() {
        float[] Cb=new float[total_nodes];
     
        for (int i=0; i<total_nodes;i++) { Cb[i]=0;}
         int[] d=new int[total_nodes];
           float[] sp=new float[total_nodes];        
        for (int s=0; s<total_nodes;s++) {           
            ArrayList<Integer> Q=new ArrayList<Integer>();
            Stack<Integer> S=new Stack<Integer>();
            ArrayList<Integer>[] P=new ArrayList[total_nodes];
            for (int i=0; i<total_nodes;i++) {
                sp[i]=0;             
                d[i]=-1;     
                P[i]=new ArrayList<Integer>();
            }
            sp[s]=1;
            d[s]=0;
            Q.add(s);
            while (!Q.isEmpty()) {
                int v=Q.remove(0); //--REmove last (dequeue)
                S.push(v);               
                for (int w:this.get_adj(v).keySet()) {                    
                       
                        if (d[w]<0) {
                            Q.add(Q.size(),w);
                            d[w]=d[v]+1;
                            
                        }
                        if (d[w]==(d[v]+1)){
                            sp[w]=sp[w]+sp[v];
                             P[w].add(v);
                        }                    
                }
            } //--End while
            //--Compute Cb
            float[] delta=new float[total_nodes];
            for (int i=0;i<total_nodes;i++) {
                //System.out.println(i+" "+sp[i]);
                delta[i]=0;
            }
            
            while(!S.isEmpty()) {
                 
                int w=S.pop();
       
                for (int v:P[w]) {                                       
                    delta[v]=delta[v]+((sp[v]/sp[w])*(1.0f+delta[w]));
                   
                }                
                if (w!=s) Cb[w]=Cb[w]+delta[w];
            } //--End S                       
        } //--End for s
        if (!directed) for (int i=0; i<total_nodes;i++) Cb[i]/=2;
        return Cb;
   }
   
    /**
    * Sarıyüce, A. E., Kaya, K., Saule, E., & Catalyürek, U. V. (2013, October). Incremental algorithms for closeness centrality. In IEEE International Conference on BigData.
    * @return array of the closeness of each node
    */
   public float[] Closeness() {
        float[] cc=new float[total_nodes];
     
        for (int i=0; i<total_nodes;i++) { cc[i]=0;}
         int[] d=new int[total_nodes];
           float[] far=new float[total_nodes];        
        for (int s=0; s<total_nodes;s++) {           
            ArrayList<Integer> Q=new ArrayList<Integer>();            
            for (int i=0; i<total_nodes;i++) {
                d[i]=-1;     
            }
            d[s]=0;
            far[s]=0;
            
            Q.add(s);
            while (!Q.isEmpty()) {
                int v=Q.remove(0); //--REmove last (dequeue)
                     
                for (int w:this.get_adj(v).keySet()) {                    
                       
                        if (d[w]<0) {
                            Q.add(Q.size(),w);
                            d[w]=d[v]+1;
                            far[s]=far[s]+d[w];
                        }                                          
                }
            } //--End while
          cc[s]=1.0f/far[s];
        } //--End for s       
        return cc;
   }
   
   public int[][] Floyd() {
       int[][] dist=new int[total_nodes][total_nodes];
       this.epath=new int[total_nodes][total_nodes];
       for (int i=0;i<total_nodes;i++) {
           for (int j=0;j<total_nodes;j++) {
               epath[i][j]=-1;
               dist[i][j]=infinity;
               if (get_adj(i).containsKey(j)) {
                   epath[i][j]=j;
                   dist[i][j]=1;
               } 
           }
           dist[i][i]=infinity;
       }
       for (int i=0;i<total_nodes;i++)
        for (int j=0; j<total_nodes;j++)
         for (int k=0; k<total_nodes;k++) {
            
                 if (dist[i][j]==dist[i][k]+dist[k][j]) {
                   ArrayList<Integer> tmp=new ArrayList<Integer>();
                     if (epath[i][j]<total_nodes) {
                        tmp.add(epath[i][j]);
                        tmp.add(epath[i][k]);
                         epath[i][j]=epath_next;
                         
                        epath_extra.put(epath_next++, tmp);
                   } else {
                       tmp=epath_extra.get(epath[i][j]);
                       tmp.add(epath[i][k]);
                       epath_extra.put(epath_next, tmp);
                   }
                 } else 
                   if (dist[i][j]>dist[i][k]+dist[k][j]) {
                     epath[i][j]=epath[i][k];
                 
                    dist[i][j]=dist[i][k]+dist[k][j];
             
                 }
             
         }   
       return dist;
   }
   
   float in_degree(int nodeid) {
       float total=0;
       for (int i:adjlist.keySet()) {
           if (adjlist.get(i).containsKey(nodeid)) total++;
       }
       return total;
   }

    float out_degree(int nodeid) {
       float total=get_adj(nodeid).size();       
       return total;
   }
  
  float density() {
      float d=this.total_edges;     
      float total_possible=(this.total_nodes*(this.total_nodes-1))/2;
      return d/total_possible;
   }
    
}

