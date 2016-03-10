package COMPONENT_GRAPHER;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
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


/**
 * Module to analyze multiple network
 * @author Etienne Lord, Jananan Pathmanathan
 * @since Mars 2016
 */
public class analyze_multiple {
    
    ////////////////////////////////////////////////////////////////////////////
    /// VARIABLES
    HashMap<Integer,HashMap<Integer,Integer>> count_edges=new HashMap<Integer,HashMap<Integer,Integer>>();
    HashMap<String,Integer> node_to_id=new HashMap<String,Integer>();
    HashMap<Integer,String> id_to_node=new  HashMap<Integer,String>(); 
    int total_nodes=0;
    int total_edges=0;
    int total_complete_edges=0;
    public int mincount=0; //minimum count to include edge in the results
    public float maxiter=0;
    util complete_network_file=new util();
    HashMap<String,Integer> total_complete_nodes=new HashMap<String,Integer>();
            
    
    ////////////////////////////////////////////////////////////////////////////
    /// FUNCTIONS
    public int getNode(String stri) {
       stri=stri.trim();
       Integer node=this.node_to_id.get(stri);
       if (node!=null) return node;
       node=this.node_to_id.size();
       node_to_id.put(stri,node);
       this.total_nodes++;
       this.id_to_node.put(node,stri);
       return node;
   }
   
   public void addEdge(int source, int dest) {
       HashMap<Integer,Integer> dest_nodes=this.count_edges.get(source);
       if (dest_nodes==null) dest_nodes=new HashMap<Integer,Integer>();
       Integer count=dest_nodes.get(dest);
       if (count==null) count=0;
       count++;
       dest_nodes.put(dest,count);
       this.count_edges.put(source, dest_nodes);
   }
   
   boolean load_network(String filename) {
       try {
           BufferedReader br=new BufferedReader(new FileReader(new File(filename)));
           String data="";
           while(br.ready()) {
               data=br.readLine();
               String[] stri=data.split("\t");
               if (stri.length>1&!data.startsWith("#")) {
                    int source=getNode(stri[0]);
                    int dest=getNode(stri[1]);
                    addEdge(source,dest); 
               }
           }
           br.close();
       } catch(Exception e) {
           return false;
       }
       total_nodes=node_to_id.size();
       return true;
   
   }
   
   public void write_network(String filename, int mincount, int type) {
       util u=new util();
       
       u.open(filename);
       u.println("#src_id\tdest_id\tedge_type\tpercent");
       for (Integer src_id:this.count_edges.keySet()) {
           HashMap<Integer,Integer> dest_ids=count_edges.get(src_id);
           for (Integer dest_id:dest_ids.keySet()) {
               int count=dest_ids.get(dest_id);
               
               if (count>=mincount) {
                   u.println(id_to_node.get(src_id)+"\t"+id_to_node.get(dest_id)+"\t"+type+"\t"+(count/maxiter));
                   
                   if (type!=4) {
                       total_complete_nodes.put(id_to_node.get(src_id),0);
                       total_complete_nodes.put(id_to_node.get(dest_id),0);
                       complete_network_file.println(id_to_node.get(src_id)+"\t"+id_to_node.get(dest_id)+"\t"+type+"\t"+(count/maxiter));
                       this.total_complete_edges++;
                   }
                   this.total_edges++;                   
               }
           }
       }
       u.close();
   }
   
    void reset() {
        this.total_nodes=0;
        this.total_edges=0;
        this.id_to_node.clear();
        this.node_to_id.clear();
        this.count_edges=new HashMap<Integer,HashMap<Integer,Integer>>();
    }
    
    void process_path(String filename, float maxiter) {
        File fi=new File(filename);
        this.maxiter=maxiter;
        complete_network_file.open(filename+"_global_complete.txt");
        complete_network_file.println("#src_id\tdest_id\tedge_type\tpercent"); 
        //        String path=fi.getAbsolutePath();
        //        filename=fi.getName();
        //        path=path.substring(0,path.length()-filename.length());
        //        
        //System.out.println(path+"\t"+filename);
        //ArrayList<String> filenames=util.listDirWithFullPath(path);
        //filenames=util.filter(filenames,filename+"_[0-9]{1,}_[1-4].txt");
        //System.out.println(filenames.size());
        // 1. For each group of graph i.e. _1, _2, _3, _4, _complete
        //    process, the network and output a new graph with all the link, their 
        //    weights, etc.
         String[] types={"_1.txt","_2.txt","_3.txt","_4.txt"};
         int[] types_s={1,2,3,4};
         for (int i=0; i<types.length;i++) {
             ArrayList<String> current_filenames=new ArrayList<String>();
             for (int j=1; j<=maxiter;j++) current_filenames.add(filename+"_"+j+types[i]);
            System.out.print("Saving type "+types_s[i]+"   - summary of "+maxiter+" variations: ");
            for (String f:current_filenames) {
               load_network(f);
           }
             write_network(filename+"_global"+types[i],this.mincount, types_s[i]);
             System.out.println(this.total_nodes+" nodes - "+this.total_edges+" edges");
             System.out.println("(saving in "+filename+"_global"+types[i]+")");
             reset();
         }
         System.out.println("Saving complete - summary of "+this.maxiter+" variations: "+this.total_complete_nodes.size()+" nodes - "+total_complete_edges+" edges");
           System.out.println("(saving in "+filename+"_global_complete.txt)");
         complete_network_file.close();
    }
    
    
    
}
