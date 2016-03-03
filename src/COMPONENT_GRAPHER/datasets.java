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


import static COMPONENT_GRAPHER.util.isNumber;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import umontreal.iro.lecuyer.util.BitVector;

/**
 * This is a new version 
 * The code to compute the graph is included in this class
 * The edge are now 2 array of integer : src_edge, dst_edge
 * @author Etienne Lord, Jananan Pathmanathan
 * @since October/November 2015, 2016
 */
public class datasets {
    
   /////////////////////////////////////////////////////////////////////////////
   /// OPTIONS 
   public double min_rand_index=0.0;
   public double min_taxa_percent=0;
   public double max_taxa_percent=1;
   public double min_taxa=1;
   
   public boolean save_graphml=false;
   public boolean nooutput=false;
   public boolean remove_undefined_column=false;
   public boolean remove_multiple_column=false;
   public boolean bipartite=false;
   public boolean save_summary=false;
   public int bipartite_type=0;
   public int bipartite_index=0;
   public int current_bipartite=0;
   public String user_state_string="";
   public HashMap<String,Integer> bipartite_node_id=new HashMap<String, Integer>();
   public String taxa="";
   
   public String tmpfile="";
   public int maxiter=1;
   util output_biparition_complete=new util();
   util output_biparition_1=new util();
   util output_biparition_2=new util();
   util output_biparition_3=new util();
   
   
   /////////////////////////////////////////////////////////////////////////////
   /// VARIABLES - Datasets
   public ArrayList<String> charlabels=new ArrayList<String>();
   public ArrayList<ArrayList<String>> statelabels=new ArrayList<ArrayList<String>>();
   public HashMap<String, Integer> index_label=new HashMap<String, Integer>();
   public ArrayList<String> label=new ArrayList<String>(); //--Taxa name
   public ArrayList<String> state=new ArrayList<String>(); //--line of char. for the taxa 
   public static Pattern isNumbers=Pattern.compile("([0-9]{1,})\\s*([0-9]{1,})");   
   public String title="";
   public String filename="";   
   public int ntax=0;
   public int nchar=0;	
   public int max_char_state=1; //--Maximum char state found for a taxon e.g. {1,2,3} =3
   public String char_matrix[][]; //--This is the original data matrixd minus the {}
   public int mode=0;
    /////////////////////////////////////////////////////////////////////////////
   /// VARIABLES - Datasets
   //--This is for rapid access to node type
   public ArrayList<HashMap<Integer,Integer>> node_id_type=new ArrayList<HashMap<Integer,Integer>>();
   
   public ArrayList<Integer> undefined_column=new ArrayList<Integer>();
   public ArrayList<Integer> multiple_column=new ArrayList<Integer>();
   
   public int char_state[]; //--char state for this partition
   /////////////////////////////////////////////////////////////////////////////
   /// STATES
   int total_state_id=0; 
   float total_states=1;
   public boolean save_inter_result=true;
   PrintWriter pw_output;
   
   ArrayList<state> states=new ArrayList<state>(); //--
   public static ArrayList<String> state_strings=new ArrayList<String>(); //--StateString for this char matrix        
   public String current_state_matrix[][]; 
   /////////////////////////////////////////////////////////////////////////////
   /// Constant
     public static String[] m_type={"","perfect concomittant","inclusion","partial concomittant","disjoint","inclusion"};
   
   /////////////////////////////////////////////////////////////////////////////
   /// Results for this datasets
 
   public static ArrayList<ArrayList<Integer>>precomp_partitions[]; //Precalculated partition   
   //public static ConcurrentHashMap <ArrayList<Integer>,Integer> vertex_numbering=new ConcurrentHashMap <>();
   //public static ConcurrentHashMap <Integer,ArrayList<Integer>> vertex_numbering_inv=new ConcurrentHashMap <>();
   //public static HashMap<edge,Integer> edges=new HashMap<>();
   ///-Edge and nodes
   //ArrayList<edge> edges=new ArrayList<>();
    // Edges are now replaced with an array of max size 
   // total_node*(total_node-1)/2
   int total_edge=0; //--To remove
   int total_edges=0;
   int[] src_edge;
   int[] dest_edge;
   int[] type_edge;
   int[] taxa_edge;
   int[] count_edge; // For multiple_edge, the stengh of the link 
   //int[] rand_edge; 
   
   ArrayList<node> nodes=new ArrayList<node>();
   
  
   //HashMap<Integer,Integer> edges_found=new HashMap<>();
   
   //--Statistics
   //HashMap<edge, Integer> edges_count=new HashMap<edge, Integer>();
   public static ConcurrentHashMap<String,Integer> identification=new ConcurrentHashMap<String,Integer>();
   public static ConcurrentHashMap<Integer,String> inv_identification=new ConcurrentHashMap <Integer,String>();
   
   int total_type0=0;
   int total_type1=0;
   int total_type2=0;
   int total_type3=0;

/////////////////////////////////////////////////////////////////////////////
   /// FUNCTIONS    
    String extract_charlabels(String line) {
        //CHARLABELS
        //		 [1] 'GEN skull, telescoping, presence'
	String num="";
	String label="";
	ArrayList<String> tmp=new ArrayList<String>();
	boolean in_num=false;
	boolean in_label=false;
	
	for (int i=0; i<line.length();i++) {
		if (line.charAt(i)=='[') {
		in_num=true;
		} else if (line.charAt(i)==']') {
		 in_num=false;
		 in_label=true;
		} else if (in_num) {
		   num+=line.charAt(i);
		} else if (in_label) {
			if (line.charAt(i)!='\'') label+=line.charAt(i);
		}
	}	
	return (label);
    }


    ArrayList<String> extract_matrix(String line) {

            String name="";
            String matrix="";
            ArrayList<String> tmp=new ArrayList<String>();

            boolean flag_in_name=false;
            // NULL CASE
            if (line.length()==0) {
                    tmp.add("");
                    tmp.add("");
                    return tmp;
            }
            // CASE WITHOUT '
            if (line.charAt(0)!='\'') {
                    flag_in_name=true;
                    for (int i=0; i<line.length();i++) {
                            if (line.charAt(i)==' '||line.charAt(i)=='\t') {
                                    flag_in_name=false;
                            } else if (flag_in_name) {
                                    name+=line.charAt(i);
                            } else if (!flag_in_name&&line.charAt(i)!=' '&&line.charAt(i)!='\t'&&line.charAt(i)!='\'') matrix+=line.charAt(i);		
                    }
            } else 
            // CASE WITH '
            if (line.charAt(0)=='\'')
                    for (int i=0; i<line.length();i++) {
                            if (line.charAt(i)=='\'') {
                                    flag_in_name=!flag_in_name;
                            } else if (flag_in_name) {
                                    name+=line.charAt(i);
                            } else if (!flag_in_name&&line.charAt(i)!=' '&&line.charAt(i)!='\t'&&line.charAt(i)!='\'') matrix+=line.charAt(i);		
                    }
            tmp.add(name);
            tmp.add(matrix);
            return tmp;		
    }
    
    int intmaxrow() {
 	return (this.state.size());
 }
    
String[][] charmatrix() {
  String  mat[][]=new String[intmaxrow()][intmaxcol()];
   
  for (int i=0; i<state.size();i++) {
	int l=0;
	boolean inside=false;
	String st="";
        
	for (int j=0;j<state.get(i).length();j++) {
		
		Character c=state.get(i).charAt(j);		
		if (inside&&(c=='}'||c==')')) {
			inside=false;
			if (st.length()>this.max_char_state) this.max_char_state=st.length();
                        mat[i][l++]=st;
			st="";
		} else 
		if (inside&&(c!=','&&c!=' ')) {
			st+=c;
		} else 
		if (c=='{'||c=='(') {
			inside=true;
			st="";
		} else
		if (!inside) {
			st+=c;
			mat[i][l++]=st;
                        st="";
		}
	}
    }
   
   return mat;
}
  // return the new matrix
  
 public Integer get_value(String s) {     
     s=s.substring(0, s.length()-1); 
     return Integer.valueOf(s.split("=")[1]);
 }
  
  public Integer get_value(String s, String id) {     
      Pattern p=Pattern.compile(id, Pattern.CASE_INSENSITIVE);
      Matcher m=p.matcher(s);
      if (m.find()) {          
          return Integer.valueOf(m.group(1));
      }
      //s=s.substring(0, s.length()-1); 
     return 0;
 }
  
     
  boolean load_morphobank_nexus(String filename) {
        ArrayList<String> tmp_state=new ArrayList<String>();
         this.filename=filename;
	 // flags!
	 boolean flag_in_matrix=false;
	 boolean flag_in_CHARLABELS=false;
	 boolean flag_in_STATELABELS=false;
        int count_section=0;
  	ArrayList<String> data=loadStrings(filename);
        if (data.isEmpty()) return false;
        for (int i=0; i<data.size();i++) {
                
                String line=data.get(i);   
                line=line.replaceAll("\t", "");                
                line=line.trim();
                //System.out.println(flag_in_matrix+" "+line);
                if (!line.isEmpty()) {
                    if ((line.indexOf("DIMENSIONS")>-1||line.indexOf("dimensions")>-1)&&line.indexOf("NTAX")>-1||line.indexOf("ntax")>-1) {
                        this.ntax=get_value(line, "ntax=([0-9]*)");
                    }
                    if (line.indexOf("BEGIN CHARACTERS;")>-1) count_section++;
                    if ((line.indexOf("DIMENSIONS")>-1||line.indexOf("dimensions")>-1)&&line.indexOf("NCHAR")>-1||line.indexOf("nchar")>-1) this.nchar=get_value(line, "nchar=([0-9]*)");;    		  
                    if (flag_in_matrix&&line.indexOf(";")>-1) flag_in_matrix=false;	
                    if (flag_in_CHARLABELS&&(line.charAt(0)==';'||line.charAt(line.length()-1)==';')) flag_in_CHARLABELS=false;	
                    if (flag_in_STATELABELS&&(line.charAt(0)==';'||line.charAt(line.length()-1)==';')) {
                            if (tmp_state.size()>1) {
                                tmp_state.remove(0);
                                statelabels.add(tmp_state);                                                
                        }
                        flag_in_STATELABELS=false;         		 	
                    }
     		// test flag	
			if (line.indexOf("CHARLABELS")>-1||line.indexOf("charlabels")>-1) { 
				flag_in_CHARLABELS=true;			
			} else 
			if (line.indexOf("STATELABELS")>-1||line.indexOf("statelabels")>-1) { 
				flag_in_STATELABELS=true;			
			} else 
			if (line.indexOf("MATRIX")>-1) { 
                            //Note: Matrix must be in Uppercase
                            flag_in_matrix=true;			
			} else 
			if (flag_in_matrix) {
				 //Try to extract
				ArrayList<String> d=extract_matrix(line);
				// Test if we have the label
                                Integer index=index_label.get(d.get(0));
                                if (d.get(1)==null||d.get(1).isEmpty()) {
                                    // Add to last state
                                    index=state.size()-1;
                                    state.set(index,state.get(index)+d.get(0));
                                         
                                } else {
                                    //Integer index=index_label.get(d.get(0));
                                    if (index==null) {
                                        index=label.size();
                                        index_label.put(d.get(0),index);
                                        label.add(d.get(0));
                                        state.add("");
                                    }
                                    state.set(index,state.get(index)+d.get(1));	
                                }
                                //
                                
                                //System.out.println(d.get(0)+"|" +index+"|"+d.get(1)+"|"+d.size()+"|"+line);
                                //System.out.println(d.get(0));
				                                  
			 }
			 else 
			if (flag_in_CHARLABELS) {
				 charlabels.add(extract_charlabels(line));
                                 if (line.indexOf('/')>0) {
                                     //--We might have state information embedded
                                 }
			 }
			else 
			if (flag_in_STATELABELS) {			
                                if (line.charAt(0)==','||line.charAt(0)==';'||line.endsWith(",")) {
				   //System.out.println(line);
                                   if (tmp_state.size()>1) {
                                       if (line.endsWith(",")) tmp_state.add(line.replaceAll("'", "").replaceAll(",",""));
                                       tmp_state.remove(0);                                   
				 	statelabels.add(tmp_state);
				 	//new state
				 	tmp_state=new ArrayList<String>();
                                   }
				 } else {
				 	tmp_state.add(line.replaceAll("'", ""));                                        
				 }
				
			 } //--End statelabels
                } //--End line empty
	} //--End data 
        //ntax=label.size();
        //nchar=intmaxcol();
        if (count_section>1) {
            System.err.println("Warning, more than one matrix found in file: "+filename);
            return false;
        }
  this.char_matrix=charmatrix();
    //--Create the different states
  try {
    if (this.nchar!=0) create_states(); 
  } catch(Exception e) {
    //e.printStackTrace();
      System.out.println("Unable to create char. state matrix.->"+this.filename);
      for (String k:this.index_label.keySet()) {
          Integer v=this.index_label.get(k);
          System.out.println(this.label.get(v)+"|"+this.state.get(v));
      }
      //printCharMatrix();
      return false;
  }
  return true;
}
  

   boolean load_simple(String filename) {  
	this.filename=filename;
       
        ArrayList<String> data=loadStrings(filename);
        if (data.isEmpty()) return false;
        for (int i=0; i<data.size();i++) {
           String line=data.get(i).trim();
           //--Do we have a phylip file ? 
           
            Matcher m=isNumbers.matcher(line);
            if (m.find()&&i==0)  {
                //--Skip line
            } else {                
                  ArrayList<String> d=extract_matrix(line);            
                  label.add(d.get(0));	
                  state.add(d.get(1));	
             }   
        }
            
        ntax=label.size();
        nchar=intmaxcol();
        this.char_matrix=charmatrix();
        
          try {
                if (this.nchar!=0) create_states(); 
              } catch(Exception e) {
                  return false;
              }
        System.out.println("Input                                : "+filename);       
        return true;   	
     }
     
  
  int intmaxcol() {   
    int c=0;
    for (int i=0; i<state.size();i++) {
            int l=0;
            boolean inside=false;
            String s=state.get(i);
            for (int j=0;j<s.length();j++) {
                    Character cs=s.charAt(j);
                    if (!inside) l++;
                    if (inside&&(cs=='}'||cs==')')) inside=false;
                    if (cs=='{'||cs=='(')  inside=true;                    
            }
            //--debug System.out.println(l+"|"+s);
            if (l>c) c=l; 
    }
    return (c);
    }

  
  public static ArrayList<String> loadStrings(String filename) {
		ArrayList<String> tmp=new ArrayList<String>();
		try {
			//Change to read UTF-8 here
			BufferedReader br=new BufferedReader(new InputStreamReader(new FileInputStream(new File(filename)),"ISO-8859-1"));
			while (br.ready()) {
                            tmp.add(br.readLine());
                        }
                        br.close();
		} catch(Exception e) {}
		return tmp;
	}  
  
        //--Here, pos is the char position
        ArrayList<String> extract_char(int pos) {
                ArrayList<String> tmp=new ArrayList<String>();
               int l=this.ntax;
               for (int i=0; i<l;i++) {		
                       tmp.add(current_state_matrix[i][pos]);
               }
               return (tmp);
       }

          //Here pos is the tax position  
          ArrayList<String> extract_char_taxa(int pos) {
                ArrayList<String> tmp=new ArrayList<String>();
               int l=this.nchar;
               for (int i=0; i<l;i++) {		
                       tmp.add(current_state_matrix[pos][i]);
               }
               return (tmp);
       }
        
      public static void CleanMemory() {
         Runtime r = Runtime.getRuntime();
         r.gc();
    }
        
       String extract_tax(int tax_pos) {
               String tmp="";
               for (int i=0;i<intmaxcol();i++) {		
                       tmp+=char_matrix[tax_pos][i];
               }
               return (tmp);
       }
 
   public static String PrintMemory() {
        String stri="System allocated memory: "+Runtime.getRuntime().totalMemory()/(1024*1024)+" MB System free memory: "+Runtime.getRuntime().freeMemory()/(1024*1024)+" MB\n"+
                    "System total core: "+Runtime.getRuntime().availableProcessors()+"\n";
         
         return stri;
    }
  
   
   public void get_info() {
       
       //--Number of column with undefined states 
       int total_undefined=0; // ? or -
       int total_multiple=0; // {1,2}
       int total_undefined_column=0;
       int total_multiple_column=0; // {1,2}
       
       for (int i=0; i<this.nchar;i++) {
           boolean found_undefined=false; // ? or -           
           boolean found_multiple=false;  // {1,2}
           for (int j=0; j<this.ntax;j++) {		
                         String s=current_state_matrix[j][i];
                         if (s.equals("?")||s.equals("-")||s.equals("*")) {
                             found_undefined=true;
                             total_undefined++;
                         }
                         if (s.length()>1) {
                             found_multiple=true;
                             total_multiple++;
                         }
           }
           if (found_undefined) {
               total_undefined_column++;
               undefined_column.add(i);
           }
           if (found_multiple) {
               total_multiple_column++;
               multiple_column.add(i);
           }
       }
       //--Output to screen some information about the char matrix       
       System.out.println("N taxa                               : "+this.ntax+" (rows)");
       System.out.println("N characters                         : "+this.nchar+" (columns)");
        //--Create the various state matrix         
       System.out.println("Total number of multistate characters: "+ this.states.size());
       System.out.println("Total number of possible variations  : "+((int)this.total_states));               
       System.out.println("Total undefined column(s)            : "+total_undefined_column);
       System.out.println("Total multiple column(s)             : "+total_multiple_column);
       System.out.println("Total undefined char                 : "+total_undefined);
       System.out.println("Total multiple char                  : "+total_multiple);
   }
   /**
    * This is the main computing routine 
    * 
    * @param mode (either 0 - taxa mode (default) or 1 (char mode)
    * @param state
    * @return 
    */
  public void compute() {       
      get_info();
      //--Clear
        for (int i=0; i<4;i++) {
           node_id_type.add(new HashMap<Integer,Integer>());
       }
      state_strings.clear();
      nodes.clear();
      identification.clear();
      inv_identification.clear();      
     
      //(min. rand index: "+this.min_rand_index+", min. shared taxa (%): "+this.min_taxa_percent);       
      System.out.println("Remove multiple state columns        : "+this.remove_multiple_column);
      System.out.println("Remove undefined columns             : "+this.remove_undefined_column);
      if (this.remove_multiple_column) this.maxiter=1;
      if (this.maxiter>1&&remove_multiple_column) System.out.println("Max. iteration (if multiple states): "+this.maxiter);      
      compute_partition();            
  }
  
  public void display_result(String filename) {
     StringBuilder st=new StringBuilder();
     st.append("===============================================================================\n");
      st.append("Results:\n");
      st.append("===============================================================================\n");
      st.append("Edges (total)                     : "+total_type0+"\n");
      int total_1=0;
      int total_2=0;
      int total_3=0;
      int unassigned_node=0;
      for (node n:nodes) if (!node_id_type.get(0).containsKey(n.id)) unassigned_node++;
      st.append("Edges type 1 (perfect)            : "+total_type1+"\n");
      st.append("Edges type 2 (inclusion)          : "+total_type2+"\n");
      st.append("Edges type 3 (partial concomitant): "+total_type3+"\n");
      st.append("\n");
      st.append("Total nodes evaluated             : "+nodes.size()+"\n");      
      st.append("Total nodes                       : "+node_id_type.get(0).size()+"\n");      
      st.append("Node (unassigned)                 : "+unassigned_node+"\n");  
      st.append("Node type 1 (perfect)             : "+node_id_type.get(1).size()+"\n");
      st.append("Node type 2 (inclusion)           : "+node_id_type.get(2).size()+"\n");
      st.append("Node type 3 (partial concomitant) : "+node_id_type.get(3).size()+"\n");
      st.append("===============================================================================\n");
      
      if (total_states>1) {
          String sti=state_strings.get(state_strings.size()-1);
          st.append("states evaluated: "+sti+"\n");
          st.append("Taxa->Character(column)|Value\n");
          st.append("----------------------------------------\n");
          
          for (int i=0; i<this.states.size();i++) {
                 state s=states.get(i);
                 //--This might fail if there is no label
                 if (this.charlabels.size()>0) {
                    st.append(this.label.get(s.pos_i)+"->"+this.charlabels.get(s.pos_j)+"|"+sti.charAt(i)+"\n");
                 } else {
                     st.append(this.label.get(s.pos_i)+"->"+(s.pos_j+1)+"|"+sti.charAt(i)+"\n");
                 }
          }
      st.append("===============================================================================\n");
   
      }
      
      System.out.println(st);     
      
      try {
          PrintWriter pw=new PrintWriter(new FileWriter(new File(filename+"_stat.txt")));
          pw.println(st);
          pw.close();
      } catch(Exception e) {}
      
      if (save_summary) {
          ArrayList<ArrayList<Integer>> CC;
          Integer[] CC_info_type1;
          boolean[] local_articulation_point=new boolean[nodes.size()];
          boolean[] global_articulation_point=new boolean[nodes.size()];
          boolean[] local_articulation_point_complete=new boolean[nodes.size()];
          boolean[] global_articulation_point_complete=new boolean[nodes.size()];
          //ArrayList<Integer>[] global_articulation_point=new ArrayList[4];
          Float[][] Triplets=new Float[4][nodes.size()];
          Float[] betweenness=new Float[nodes.size()];
          Float[] closeness=new Float[nodes.size()];
          Float[] in_degree2=new Float[nodes.size()];
          Float[] in_degree2_norm=new Float[nodes.size()];
          Float[] path_len4_type3=new Float[nodes.size()];
          Float[] path_loop_len4_type3=new Float[nodes.size()];
          int[][] degrees=new int[this.nodes.size()][6];
          CC_info_type1=new Integer[this.nodes.size()];
          
          Integer[] CC_info_complete=new Integer[this.nodes.size()];
          ArrayList<String>[] Progressive_transition=new  ArrayList[this.nodes.size()];
          Integer[] max_sp_type3=new Integer[this.nodes.size()];
          Integer[] max_sp_complete=new Integer[this.nodes.size()];
          //CC_info=new Integer[this.nodes.size()];
          //int CC_node_complete=0;
          float total_triplet_type3=0;
           float total_triplet_complete=0;
           float total_triplet_type2=0;
           int total_CC_type1=0;
           int total_CC_complete=0;
           
          for (int i=0; i<this.nodes.size();i++) {
              for (int j=0; j<6;j++) degrees[i][j]=0;
              CC_info_type1[i]=0;
              CC_info_complete[i]=0;
              max_sp_type3[i]=0;
              max_sp_complete[i]=0;
              in_degree2_norm[i]=0.0f;
              in_degree2[i]=0.0f;
              path_len4_type3[i]=0.0f;
              path_loop_len4_type3[i]=0.0f;
              Progressive_transition[i]=new ArrayList<String>(0);
          }
//              local_articulation_point[i]=new ArrayList<Integer>(); //--Node which are articulation point for this type
//              global_articulation_point[i]=new ArrayList<Integer>(); //--Node which are articulation point for this type
//              //Triplets[i]=new ArrayList<Float>(); //--Node which are articulation point for this type
//              //betweenness[i]=new ArrayList<Float>();
//          }
          
          //--Analyse graph
          
         
           System.out.println("================================ SUMMARY ======================================");
           System.out.println("Network\tVertex\tEdges\tDensity");
           System.out.println("----------------------------------------");
          
           for (int type=0; type<4;type++) {
              graph g=new graph();
              //--Add edge
              for (int i=0; i<this.total_edge;i++) {
                  if (type==0) {
                      int src=g.addNode(this.src_edge[i]);
                      int dest=g.addNode(this.dest_edge[i]);
                      switch(this.type_edge[i]) {
                          case 1: degrees[this.src_edge[i]][0]++;
                                  degrees[this.dest_edge[i]][0]++;
                                  degrees[this.src_edge[i]][1]++;
                                  degrees[this.dest_edge[i]][1]++;                                  
                                  break;                                                        
                          case 2: degrees[this.src_edge[i]][3]++;                                                                    
                                  degrees[this.dest_edge[i]][2]++;                                  
                                  break;                                  
                          case 3: degrees[this.src_edge[i]][4]++;
                                  degrees[this.dest_edge[i]][4]++;
                                  degrees[this.src_edge[i]][5]++;
                                  degrees[this.dest_edge[i]][5]++;                                  
                                  break;                                 
                      }
                      g.addEdge(src, dest);
                      g.addEdge(dest, src);                      
                      g.total_edges++;
                      g.directed=false;                      
                      
                  }
                  if (type==1&&this.type_edge[i]==1) {
                      int src=g.addNode(this.src_edge[i]);
                      int dest=g.addNode(this.dest_edge[i]);
                      g.addEdge(src, dest);
                      g.addEdge(dest, src);                      
                        g.total_edges++;
                      g.directed=false;            
                     
                  }
                  if (type==2&&this.type_edge[i]==2) {
                      int src=g.addNode(this.src_edge[i]);
                      int dest=g.addNode(this.dest_edge[i]);
                      g.addEdge(src, dest);
                      g.total_edges++;
                      g.directed=true;
                  } 
                  if (type==3&&this.type_edge[i]==3) {
                      int src=g.addNode(this.src_edge[i]);
                      int dest=g.addNode(this.dest_edge[i]);
                      g.addEdge(src, dest);
                      g.addEdge(dest, src);   
                      g.total_edges++;
                      g.directed=false;                        
                  } 
              }
              g.total_nodes=g.id_to_old_id.size();
                            
              System.out.println(""+(type==0?"complet":type)+"\t"+g.total_nodes+"\t"+g.total_edges+" \t"+g.density());
              //--calculate some statistic on each graph               
             
              if (type==0) {
                   for (int i=0; i<g.total_nodes;i++) {                          
                     Triplets[0][g.id_to_old_id.get(i)]=g.find_triplet(i, false);                     
                    total_triplet_complete+= Triplets[0][g.id_to_old_id.get(i)];
                    if (Triplets[0][g.id_to_old_id.get(i)]==0)Triplets[0][g.id_to_old_id.get(i)]=null; //--for display purpose
                     local_articulation_point_complete[g.id_to_old_id.get(i)]=g.is_local_articulation_point(i);
                     global_articulation_point_complete[g.id_to_old_id.get(i)]=g.is_global_articulation_point(i);
                   }
                   ArrayList<ArrayList<Integer>>tmp=g.getCC();
                    total_CC_complete=tmp.size();
                    for (int i=0; i<tmp.size();i++) {
                      ArrayList<Integer>cc=tmp.get(i);
                      for (int w:cc) CC_info_complete[g.id_to_old_id.get(w)]=i+1;
                  }
                   
              }
              if (type==1) {
                  CC=g.getCC();
                  total_CC_type1=CC.size();
                  //--REname each node in cc
                  for (int i=0; i<CC.size();i++) {
                      ArrayList<Integer>cc=CC.get(i);
                      for (int w:cc) CC_info_type1[g.id_to_old_id.get(w)]=i+1;
                  }
              }              
              if (type==2) {
                 for (int i=0; i<g.total_nodes;i++) {                          
                   
                   in_degree2[g.id_to_old_id.get(i)]=g.in_degree(i);
                   in_degree2_norm[g.id_to_old_id.get(i)]=  in_degree2[g.id_to_old_id.get(i)]/(float)(g.total_nodes-1);
                }   
              }
              if (type==3) {
                  float[] tmp=g.Betweenness();
                  float[] tmp_c=g.Closeness();
                  for (int i=0; i<g.total_nodes;i++) {
                  betweenness[g.id_to_old_id.get(i)]=tmp[i];
                  closeness[g.id_to_old_id.get(i)]=tmp_c[i];
                  local_articulation_point[g.id_to_old_id.get(i)]=g.is_local_articulation_point(i);
                  global_articulation_point[g.id_to_old_id.get(i)]=g.is_global_articulation_point(i);
                  graph.results r=g.findLoops(i,false);
                  path_len4_type3[g.id_to_old_id.get(i)]=r.total_len4;
                  path_loop_len4_type3[g.id_to_old_id.get(i)]=r.total_loop3+r.total_loop4;
                  
                    //if (g.is_local_articulation_point(i)) local_articulation_point[0].add(g.id_to_old_id.get(i));
                  //if (g.is_global_articulation_point(i)) global_articulation_point[0].add(g.id_to_old_id.get(i)); 
                  Triplets[3][g.id_to_old_id.get(i)]=g.find_triplet(i, false);
                  total_triplet_type3+= Triplets[3][g.id_to_old_id.get(i)];
                }
              }
          }
          
          //--Test for transitive progress (Shortpath>2 in type 3 but never smaller in complete
          //--Note: the graph is mixed
          graph g3=new graph();
          graph gc=new graph();
          for (int i=0; i<this.total_edge;i++) {
                      int type=this.type_edge[i];
                      int src=gc.addNode(this.src_edge[i]);
                      int dest=gc.addNode(this.dest_edge[i]);                      
                      gc.addEdge(src, dest);
                      if (type!=2)gc.addEdge(dest, src);                                              
                      src=g3.addNode(this.src_edge[i]);
                      dest=g3.addNode(this.dest_edge[i]); 
                      if (type==3) {
                          g3.addEdge(dest, src);
                          g3.addEdge(src, dest);
                      }
                      
                     
         }
          g3.directed=false; 
          gc.directed=true;
          g3.total_nodes=g3.id_to_old_id.size();
          gc.total_nodes=gc.id_to_old_id.size();
          //--Now execute Floyd
          int[][] gcf= gc.Floyd();
          int[][] g3f=g3.Floyd();
          //--Search for progressive_tr
          for (int i=0; i<g3.total_nodes;i++) {
              int original_id=g3.id_to_old_id.get(i);
              max_sp_type3[original_id]=-1; 
              max_sp_complete[original_id]=-1;
              int original_id_gc=gc.old_id_to_id.get(original_id);
              for (int j=0; j<g3.total_nodes;j++) {
              if (i!=j) {
                int original_id_j=g3.id_to_old_id.get(j);
                //int nodeid_g3=i;
                int nodeid_gc=gc.old_id_to_id.get(original_id_j);
                Integer len_g3=g3f[i][j];
                Integer len_gc=gcf[original_id_gc][nodeid_gc];
                    //System.out.println(original_id+" "+original_id_j+" "+len_g3+" "+len_gc);
                    if (len_g3<graph.infinity) {
                        if (len_g3>2&&len_gc>=len_g3) {
                            //System.out.println(nodes.get(original_id).complete_name+" -> "+nodes.get(original_id_j).complete_name+ " "+len_g3+" "+len_gc);
                            Progressive_transition[original_id].add(nodes.get(original_id_j).complete_name);
                        }
                    }
                   if (len_g3>max_sp_type3[original_id]&&len_g3<graph.infinity) max_sp_type3[original_id]=len_g3;
                   if (len_gc>max_sp_complete[original_id]&&len_gc<graph.infinity) max_sp_complete[original_id]=len_gc;
              }
          }
         }
          //--Do we test for taxa in node noame
          ArrayList<Integer> taxa_pos=new ArrayList<Integer>();
          boolean number=false;
          if (this.taxa.length()!=0) {              
               String[] t=this.taxa.split(",");
               try {
                   for (String s:t) {   
                      Integer i=Integer.valueOf(s);
                      taxa_pos.add(i); //--Add position starting at 1
                   }
                   number=true;
                 } catch(Exception e){}
              if (!number) {
                  //--Try to match taxa position
                  for (String s:t) {
                      for (int i=0; i<this.label.size();i++) {
                          if (this.label.get(i).toLowerCase().indexOf(s)>-1) taxa_pos.add(i+1);
                      }
                  }
                  
              }
          }
          
          //--output node statistics in/out degress
          util u4=new util();
          System.out.println("----------------------------------------");
          System.out.println("\nSaving vertex degree information to "+filename+"_degrees.txt");
          u4.open(filename+"_degrees.txt");
          u4.println("nodeid\tin_degree1\tout_degree1\tin_degree2\tout_degree2\tin_degree3\tout_degree3");
          for (node n:this.nodes) {
              u4.println(n.id+"\t"+degrees[n.id][0]+"\t"+degrees[n.id][1]+"\t"+degrees[n.id][2]+"\t"+degrees[n.id][3]+"\t"+degrees[n.id][4]+"\t"+degrees[n.id][5]);
          }
          u4.close();
          try {
              System.out.println("Saving summary information to "+filename+"_summary.txt");
              PrintWriter pw=new PrintWriter(new FileWriter(new File(filename+"_summary.txt")));
              // Output informations
              pw.println("nodeid\tContains "+taxa+"\ttype_1\ttype_2\ttype_3\ttype_complete\tchar.\tstate\tchar.|states\tCC_type1\tCC_complete\tlocal_ap_type3\tglobal_ap_type3\tlocal_ap_complete\tglobal_ap_complete\tin_degree_type2\tNorm._indegree_type2\tBetweenness_type3\tCloseness_type3\tTriplet_type3\t%_Triplet_type3\tTriplet_complete\tt%_Triplet_complete\tMax_shortest_path_type3\tMax_shortest_path_complete\tConvergence\tProgressive_Transition_total\tProgressive_Transition_end_node\tContains\tPercent_contained\t"+taxa+"\tTaxa");
               //--Some counter
              int total_taxa=0;
              int total_ap_local_type3=0;
              int total_ap_local_complete=0;
              int total_ap_global_type3=0;
              int total_ap_global_complete=0;
              int total_progressive=0;
              for (node n:this.nodes) {
                 //--Get the CC for this node
                  String taxas=this.get_taxa(n.identification).toLowerCase();                                    
                  ArrayList<Integer> contains_total=(ArrayList<Integer>)util.intersection(n.identification,taxa_pos);                  
                  boolean contain_taxa=!contains_total.isEmpty();
                  float percent_contains=contains_total.size()*100/n.identification.size();
                  
                  //taxas.indexOf(this.taxa.toLowerCase())>-1);
                  //if (this.taxa.length()==0) contain_taxa=false;
                  //--Some counter
                  if (contain_taxa)total_taxa++;
                 if (local_articulation_point[n.id]) total_ap_local_type3++;
                 if (global_articulation_point[n.id]) total_ap_global_type3++;
                 if (local_articulation_point_complete[n.id]) total_ap_local_complete++;
                 if (global_articulation_point_complete[n.id]) total_ap_global_complete++;
                 
                 
                 String max_sp3="Inf.";
                 if (max_sp_type3[n.id]!=null&&max_sp_type3[n.id]<graph.infinity) max_sp3=""+max_sp_type3[n.id];
                 if (max_sp_type3[n.id]==null||max_sp_type3[n.id]<1) max_sp3="";
                 String max_spc="Inf.";
                 if (max_sp_complete[n.id]!=null&&max_sp_complete[n.id]<graph.infinity) max_spc=""+max_sp_complete[n.id];
                 if (max_sp_complete[n.id]==null||max_sp_complete[n.id]<1) max_spc="";
                 
                  total_progressive+=Progressive_transition[n.id].size();
                  pw.println(n.id+"\t"+(contain_taxa?"x":" ")+"\t"+(node_id_type.get(1).containsKey(n.id)?"x":" ")+"\t"+(node_id_type.get(2).containsKey(n.id)?"x":" ")+
                         "\t"+(node_id_type.get(3).containsKey(n.id)?"x":" ")+"\t"+(node_id_type.get(0).containsKey(n.id)?"x":" ")+"\t"+n.index+"\t"+n.state_matrix+"\t"+n.complete_name+"\t"+
                          (CC_info_type1[n.id]==null?" ":CC_info_type1[n.id])+"\t"+
                          (CC_info_complete[n.id]==null?" ":CC_info_complete[n.id])+"\t"+
                          (local_articulation_point[n.id]?"x":" ")+"\t"+(global_articulation_point[n.id]?"x":" ")+"\t"+
                          (local_articulation_point_complete[n.id]?"x":" ")+"\t"+(global_articulation_point_complete[n.id]?"x":" ")+"\t"+
                         (in_degree2[n.id]>0?in_degree2[n.id]:" ")+"\t"+(in_degree2[n.id]>0?in_degree2_norm[n.id]:" ")+"\t"+
                          (betweenness[n.id]==null?" ":betweenness[n.id])+"\t"+
                          (closeness[n.id]==null?" ":closeness[n.id])+"\t"+
                         (Triplets[3][n.id]==null?" ":Triplets[3][n.id])+"\t"+(Triplets[3][n.id]==null?" ":Triplets[3][n.id]*100.0/total_triplet_type3)+"\t"+                         
                         (Triplets[0][n.id]==null?" ":Triplets[0][n.id])+"\t"+(Triplets[0][n.id]==null?" ":Triplets[0][n.id]*100.0/total_triplet_complete)+"\t"+
                          max_sp3+"\t"+max_spc+"\t"+
                          (path_loop_len4_type3[n.id]==null||path_len4_type3[n.id]==null?" ":path_loop_len4_type3[n.id]/path_len4_type3[n.id])+"\t"+                         
                          Progressive_transition[n.id].size()+"\t"+Progressive_transition[n.id]+"\t"+
                          (contain_taxa?"x":" ")+"\t"+ percent_contains+"\t"+this.get_taxa(n.identification)
                         
                 );
              }
              //--Summary
              pw.println("Total\t"+total_taxa+"\t"+this.node_id_type.get(1).size()+"\t"+this.node_id_type.get(2).size()+"\t"+this.node_id_type.get(3).size()+"\t"+this.node_id_type.get(0).size()+"\tNA\tNA\tNA\t"+total_CC_type1+"\t"+total_CC_complete+"\t"+total_ap_local_type3+"\t"+total_ap_global_type3+"\t"+total_ap_local_complete+"\t"+total_ap_global_complete+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"+total_progressive+"\t\t"+total_taxa);
              pw.println("Note: 'x' indicates presence, CC stands for Connected Components, local_ap stands for Local Articulation Point, global_ap stands for Global Articulation Point, Triplets stands for linear series of 3 nodes where terminal nodes are not connected, Convergence stands for the ratio of loops found in paths of lenght <= 4 from the starting nodes, Progressive convergence stand for short path of length>2 in type 3 network that is not smaller in complet network.");
              pw.flush();
              pw.close();
              System.out.println("===============================================================================");
          } catch(Exception e) {e.printStackTrace();}
      }
      
  }
  
   
  //Type is the new mode October 2015
    // This include:
    // 1. the BitVector
    // 2. Computation by nodes and not by the partition
   
  public boolean compute_partition() {
     
       ///////////////////////////////
      // Precompute nodes id  
        int total_column=0;
      for (int i=0; i<this.nchar;i++) {
         
      if (!(remove_multiple_column&&multiple_column.contains(i))&&!(remove_undefined_column&&undefined_column.contains(i))) {
         total_column++;
          ArrayList<String>d=get_states_column(i);
          
          for (String s:d) { 
            int index=0;
              if (!s.equals("?")&&!s.equals("-")&&!s.equals("*")) {
                String node_name=""+(i+1)+"_"+s;
                String state=s;
                if (util.isNumber(s)) {
                    index=Integer.valueOf(s); 
                    try {
                        if (this.statelabels.size()>0&&index<this.statelabels.size()) {
                          state=this.statelabels.get(i).get(index);
                      }
                      } catch(Exception ez){
                          //ez.printStackTrace();
                         System.out.println("No information for:"+node_name+" "+i+" "+index+this.statelabels.get(i)+" "+s);
                      }
                }
                //--Create the node
                //System.out.println(node_name);
                if (!identification.contains(node_name)) {
                    node n=new node(node_name,nodes.size());
                    if (this.charlabels.size()>0) {                        
                        n.complete_name=charlabels.get(i)+"|"+state;
                        n.char_label=charlabels.get(i);                                                
                    } else {
                        n.complete_name="char. "+(i+1)+"|"+state;
                        n.char_label="char. "+(i+1);
                    }
                    n.state_label=state;                                           
                    n.state_matrix=s;     
                    n.index=(i+1);
                    n.multistate=get_node_multistate(n);
                    nodes.add(n);                    
                    identification.put(n.name,n.id);
                    inv_identification.put(n.id, n.name);
                }                
            }            
           }  
         }
      }
     
      
      System.out.println("Total treated column                 : "+total_column);
      System.out.println("Total possible nodes                 : "+nodes.size());
      System.out.println("===============================================================================");        
//--Allocate new edges
        int max_iter=this.maxiter;
      if (this.total_states==1) max_iter=1;
      if (this.total_states<1000&&this.maxiter>this.total_states) max_iter=(int)this.total_states;
        this.total_edges=(((nodes.size()*(nodes.size()-1))/2)*max_iter);
        
        System.out.println("Trying to allocate memory for "+total_edges+" possible edges.");
        src_edge=new int[total_edges];
        dest_edge=new int[total_edges];
        type_edge=new int[total_edges];
        taxa_edge=new int[total_edges];
        count_edge=new int[total_edges];
        ///////////////////////////////////
        /// Allocate edge and set as unset
       int l2=nodes.size();    
              int p=0;
              for (int i=0; i<l2;i++) {
                  for (int j=0; j<l2;j++) {
                      if (i<j) {
                        src_edge[p]=nodes.get(i).id;                              
                        dest_edge[p]=nodes.get(j).id;
                        count_edge[p]=0;
                        type_edge[p]=-1;
                        taxa_edge[p]=-1;
                        p++;
                      }
                  }
              }
       ///////////////////////////////
      // Precompute partition      
      System.out.println("Computing partition...");
      precomp_partitions=new ArrayList[this.nchar];
     
      for (int state_id=0; state_id<max_iter;state_id++) {         
          
          String solution=prepare_current_state_matrix(state_id, false);
          state_strings.add(solution);
            ///////////////////////////////
      // Bipartition
       
            String f=filename+".bipartite";
           String f2=filename+"_complete_"+solution+".txt";
                   
            
            String f_complete=f+"_complete_"+solution+".txt";
            String f_1= f+"_"+solution+"_1.txt";
            String f_2= f+"_"+solution+"_2.txt";
            String f_3= f+"_"+solution+"_3.txt";
            String f_id= f+"_"+solution+"_id.txt";          
            if (bipartite) {
                output_biparition_complete.open(f_complete);
            
                output_biparition_1.open(f_1);
                output_biparition_2.open(f_2);
                output_biparition_3.open(f_3);
            
                bipartite_index=nodes.size(); //starting bipartition id
            }
          System.out.println("===============================================================================");
          System.out.println("Current iteration : "+(state_id+1)+"/"+max_iter+ " (states: "+(solution.equals("")?"none":solution)+")\n(saving to: "+f2+")");
          System.out.println("===============================================================================");
         
         
        //--Precompute the partition
         //Reset OR not? TO DO. IF maxiter is not set, we should pile up results?
            
         this.total_edge=0; 
         this.node_id_type.get(0).clear();
         this.node_id_type.get(1).clear();
         this.node_id_type.get(2).clear();
         this.node_id_type.get(3).clear();
          for (int i=0; i<nodes.size();i++) {                       
              //--This is new, we do it for the node
              node n=nodes.get(i);                            
              ArrayList<String> stris=extract_char(n.index-1);        
              n.partition=new BitVector(stris.size());
              for (int j=0; j<stris.size();j++) {
                  if (stris.get(j).equals(n.state_matrix)) {
                      n.identification.add(j+1); //To have the taxa numbering starting at 1
                      n.partition.setBool(j, true);
                  }
              }
              n.total_taxa=util.total_bitset(n.partition);
              nodes.set(i, n);
             }                  
            int l=nodes.size();    
            int total=l*(l-1)/2;
            int total_10p=total/10;
              int k=0;
              long timerunning=System.currentTimeMillis();              
              for (int i=0; i<l;i++) {
                  for (int j=0; j<l;j++) {
                      if (i<j&&nodes.get(i).index!=nodes.get(j).index) {                                                            
                                compute_persistance_and_bipartition(nodes.get(i), nodes.get(j));
                                k++;
                                
                                if (k%total_10p==0) {
                                    long elapsed=System.currentTimeMillis()-timerunning;
                                    System.out.println(k+" / "+total+" ( "+util.msToString(elapsed)+")");                                    
                                    CleanMemory();
                                }
                      }
                  }
              }  //--End i
                long elapsed=System.currentTimeMillis()-timerunning;
               System.out.println(total+" / "+total+" ( "+util.msToString(elapsed)+")");    
             if (this.bipartite) {
               output_biparition_complete.close();
                output_biparition_1.close();
                 output_biparition_2.close();
                output_biparition_3.close();     
             }
                if (!nooutput) export_edgelist(filename+"_"+solution);   
                if (save_graphml) {           
                     export_graphml(filename+"_"+solution+"_complete",0);
                     export_graphml(filename+"_"+solution+"_1",1);
                     export_graphml(filename+"_"+solution+"_2",2);
                     export_graphml(filename+"_"+solution+"_3",3);
                }
                display_result(filename+"_"+solution);
                if (this.bipartite) { 
                    System.out.println("=============================== BIPARTITION ===================================");                
                    System.out.println("Saving bipartition files to : "+f_complete);
                    System.out.println("Saving bipartition node identification to : "+f_id);
                    output_biparition_complete.open(f_id);
                    for (String m:bipartite_node_id.keySet()) {
                        String taxa=get_taxa(m);
                        output_biparition_complete.println(bipartite_node_id.get(m)+"\t\""+taxa+"\"");
                    }
                    for (node n:nodes) {
                        output_biparition_complete.println(n.id+"\t\""+n.complete_name+"\"");
                    }
                    output_biparition_complete.close();
                    System.out.println("===============================================================================");                
                 
                } //End saving bipartite graph
                
            } //--End state
                 return true;
    }    
    
          public void prune(int min_edge) {
              for (node n:nodes) {
                  if (n.edgecount<=min_edge) {
                      node_id_type.get(0).remove(n.id);
                      node_id_type.get(1).remove(n.id);
                      node_id_type.get(2).remove(n.id);
                      node_id_type.get(3).remove(n.id);
                  }
              }
              
              
          }

//          //--Get the numbering for this vertex
//          public static int get_vertex(ArrayList<Integer> query) {
//              if (vertex_numbering.containsKey(query)) return vertex_numbering.get(query);     
//              int n=vertex_numbering.size();
//              vertex_numbering.put(query, n);
//              vertex_numbering_inv.put(n,query);
//              return n;
//          }

          public void printCharMatrix() {
              for (int i=0;i<ntax;i++) {
                  System.out.print(label.get(i)+"\t");
                  for (int j=0; j<nchar;j++) {
                      System.out.print(char_matrix[i][j]+"\t"); 
                  }
                  System.out.println("");
      }
  }
  
  public void find_duplicate_partition() {
      int len=precomp_partitions.length;
      HashMap<ArrayList<ArrayList<Integer>>,Integer> tmp=new HashMap<ArrayList<ArrayList<Integer>>,Integer>();
      
      for (int i=0; i<len;i++) {
          ArrayList<ArrayList<Integer>> p=precomp_partitions[i];
          if (tmp.containsKey(p)) {
              tmp.put(p, tmp.get(p)+1);
          } else {
              tmp.put(p, 1);
          }          
      }
      int c=0;
      for (int i:tmp.values()) {
          if (i>2) {
              System.out.println(i);
              c+=i;
          }
      }
      System.out.println(c);
     
      
  }
  
  //--Separate the graph 
   public void export_edgelist(String filename) {
       //  noeud1 separateur noeud2 separateur partition_de_taxa_commune separateur directed(or undirected) separateur type_d_arete(1,23)
       try {            
            //-Type 0
            PrintWriter pw=new PrintWriter(new FileWriter(new File(filename+"_complete.txt")));                                       
            pw.println("#src_id\tdest_id\tedge_type\tnumber_common_taxa");
            for (int i=0; i<this.total_edge;i++) {               
                    if (type_edge[i]!=-1) pw.println(""+src_edge[i]+"\t"+dest_edge[i]+"\t"+type_edge[i]+"\t"+taxa_edge[i]);                     
                }    
            pw.close();   
            //-Type 1
            pw=new PrintWriter(new FileWriter(new File(filename+"_1.txt")));                       
                pw.println("#src_id\tdest_id\tedge_type\tnumber_common_taxa"); 
                for (int i=0; i<total_edge;i++) {               
                   if (type_edge[i]==1) pw.println(""+src_edge[i]+"\t"+dest_edge[i]+"\t"+type_edge[i]+"\t"+taxa_edge[i]);                     
                }  
            pw.close();            
            //--Type 2
            pw=new PrintWriter(new FileWriter(new File(filename+"_2.txt")));                       
                pw.println("#src_id\tdest_id\tedge_type\tnumber_common_taxa");
                for (int i=0; i<total_edge;i++) {               
                  if (type_edge[i]==2) pw.println(""+src_edge[i]+"\t"+dest_edge[i]+"\t"+type_edge[i]+"\t"+taxa_edge[i]);                                        
                }  
            pw.close();                
//            //-Type 3
            pw=new PrintWriter(new FileWriter(new File(filename+"_3.txt")));                       
                pw.println("#src_id\tdest_id\tedge_type\tnumber_common_taxa");
                for (int i=0; i<total_edge;i++) {               
                   if (type_edge[i]==3) pw.println(""+src_edge[i]+"\t"+dest_edge[i]+"\t"+type_edge[i]+"\t"+taxa_edge[i]);                                 
                }  
            pw.close();     
            //--Dict
            Collections.sort(nodes);
             pw=new PrintWriter(new FileWriter(new File(filename+"_id.txt")));                                       
                for (node n:nodes) pw.println(n.id+"\t"+n.complete_name+"\t"+n.char_label+"\t"+n.state_label+"\t"+n.state_matrix+"\t"+n.edgecount+"\t"+n.in_edgecount+"\t"+n.out_edgecount+"\t"+n.identification.size()+"\t"+get_taxa(n.identification));
            pw.close();

        } catch(Exception ex) {ex.printStackTrace();}
        
    }

    @Override
    public String toString() {
        return "ntax: "+this.ntax+" nchar: "+this.nchar;
    }
    
  
   
   public void export_graphml(String filename, int type) {
       //System.out.println("**"+this.charlabels.size());   
      
       try {      
             // for (String s:this.charlabels) System.out.println(s);
             PrintWriter pw=new PrintWriter(new FileWriter(new File(filename+".graphml")));     
               pw.println("<?xml version='1.0' encoding='UTF-8' standalone='no'?>");
               pw.println("<graphml>");
               //--Attributes
               //pw.println("<key id='k1' for='edge' attr.name='weight' attr.type='double'/>");
               pw.println("<key id='k2' for='edge' attr.name='type' attr.type='double'/>");
               //pw.println("<key id='r1' for='edge' attr.name='randindex' attr.type='double'/>");
               pw.println("<key id='k1' for='edge' attr.name='total_shared_taxa' attr.type='double'/>");
               
               pw.println("<key attr.name='interaction' attr.type='string' for='egde' id='interaction'/>");             
               pw.println("<key id='k3' for='node' attr.name='fullname' attr.type='string'/>");
               pw.println("<key id='k4' for='node' attr.name='number_of_taxa' attr.type='double'/>");               
               pw.println("<key id='k5' for='node' attr.name='partition' attr.type='string'/>");
               pw.println("<key id='k6' for='node' attr.name='total_edges' attr.type='double'/>");
               pw.println("<key id='k61' for='node' attr.name='in_edges' attr.type='double'/>");
               pw.println("<key id='k62' for='node' attr.name='out_edges' attr.type='double'/>");
               pw.println("<key id='k7' for='node' attr.name='associated_character_column' attr.type='double'/>");
               pw.println("<key id='k8' for='node' attr.name='charlabel' attr.type='string'/>");
               pw.println("<key id='k9' for='node' attr.name='statelabel' attr.type='string'/>");
               pw.println("<key id='k10' for='node' attr.name='statematrix' attr.type='string'/>");
               
               pw.println("<key id='k11' for='node' attr.name='total_taxa' attr.type='double'/>");
               pw.println("<key id='k12' for='node' attr.name='taxa_id' attr.type='string'/>");
               
               //if (type==1) {
                pw.println("<graph edgedefault='undirected' id='"+this.title+"'>"); 
               //} else {
               //    pw.println("<graph edgedefault='directed' id='"+this.title+"'>"); 
               //}
               HashMap<Integer,Integer> nodes_id=node_id_type.get(type);
               for (node n:nodes) {                   
                   if (nodes_id.containsKey(n.id)) {
                    pw.println("<node id='"+n.name+"'>");
                    pw.println("<data key='k3'>"+n.complete_name+"</data>");
                    pw.println("<data key='k4'>"+n.count+"</data>");                    
                    pw.println("<data key='k6'>"+n.edgecount+"</data>");
                    pw.println("<data key='k7'>"+n.index+"</data>");
                    pw.println("<data key='k8'>"+n.char_label+"</data>");
                    pw.println("<data key='k9'>"+n.state_label+"</data>");
                    pw.println("<data key='k10'>"+n.state_matrix+"</data>");
                    pw.println("<data key='k61'>"+n.in_edgecount+"</data>");
                    pw.println("<data key='k62'>"+n.out_edgecount+"</data>");
                    pw.println("<data key='k11'>"+n.total_taxa+"</data>");
                    pw.println("</node>");
                   }
               }      
               
               for (int i=0; i<total_edge;i++) {
                  
                  if (type==type_edge[i]||type==0) {
                   if (type==1||type==3) {
                       pw.println("<edge directed='false' source='"+inv_identification.get(src_edge[i])+"' target='"+inv_identification.get(dest_edge[i])+"'>");
                   } else {
                       pw.println("<edge directed='true' source='"+inv_identification.get(src_edge[i])+"' target='"+inv_identification.get(dest_edge[i])+"'>");
                   }
                    pw.println("<data key='k1'>"+taxa_edge[i]+"</data>");                    
                    pw.println("<data key='k2'>"+type_edge[i]+"</data>");
                    //--TO DO
                    //pw.println("<data key='r1'>"+e.randindex+"</data>");
                    pw.println("<data key='interaction'>"+m_type[type_edge[i]]+"</data>");
                   pw.println("</edge>");
                  }
                  // if (e.type!=1) pw.println("<edge directed='true' source='"+e.source_str+"' target='"+e.dest_str+"'/>");
                   //if (e.type==1) pw.println("<edge  directed='false' source='"+e.source_str+"' target='"+e.dest_str+"'/>");
               }
               pw.println("</graph>");
               pw.println("</graphml>");
               pw.close();
           } catch(Exception e) {e.printStackTrace();}
   }
   

   
     /**
      * This must be called after the preparation of the char matrix
      */
     public void create_states() {
         ArrayList<String> st=new ArrayList<String>();
         int total=0;
         this.current_state_matrix=new String[this.ntax][this.nchar];        
         for (int i=0; i<this.ntax;i++) {
             for (int j=0; j<this.nchar;j++) {
                 this.current_state_matrix[i][j]=char_matrix[i][j];
                 if (char_matrix[i][j].length()>1) {
                     state s=new state();
                     s.pos_i=i;
                     s.pos_j=j;
                     s.state_id=this.total_state_id++;
                     s.state=char_matrix[i][j];                     
                     this.total_states*=s.state.length();
                     states.add(s);                    
                     st.add(s.state);
                 }
             }
         }          
     }
     
     public String prepare_current_state_matrix(int state_id, boolean rand) {
          String current_state=""; 
         boolean ok=false;                  
         if (!rand&&total_states<1000) {
             //--Generate the combinations if not generated
            if (state_strings.isEmpty()) {
                ArrayList<String> input=new ArrayList<String>();
               for (state s:states) input.add(s.state);
                state_strings=util.combinations(input);
            }
             current_state=state_strings.get(state_id);
             for(int i=0; i<states.size();i++) {
                  state s=states.get(i);
                  this.current_state_matrix[s.pos_i][s.pos_j]=""+current_state.charAt(i);
              }
         } else if (!user_state_string.isEmpty()&&user_state_string.length()>=states.size()) {
             current_state=user_state_string;             
             for(int i=0; i<states.size();i++) {
                  state s=states.get(i);
                  this.current_state_matrix[s.pos_i][s.pos_j]=""+current_state.charAt(i);
              }
         } else {
         
            while (!ok) {
                current_state="";         
                for (int i=0; i<states.size();i++) {
                      state s=states.get(i);
                      //--Randomly pick a state
                      Random r=new Random();
                      int pos=r.nextInt(s.state.length());
                      current_state+=s.state.charAt(pos);
                      this.current_state_matrix[s.pos_i][s.pos_j]=""+s.state.charAt(pos);
                  }
                if (!state_strings.contains(current_state)) {                   
                    ok=true;
                }
            } 
         }
         return current_state;
     }
     
     //--Calculate the fitness of this iteration
     public float fitness() {
         return 0.0f;
     }
          /**
      * This return the possible state found in a column
      * @param j
      * @return 
      */
     ArrayList<String> get_column(int j) {
         ArrayList<String> temp=new ArrayList<String>();
         for (int i=0; i<this.ntax;i++) {
             String d=this.char_matrix[i][j];
             temp.add(""+d);            
         }
         return temp;
     }
     
     /**
      * This return the possible state found in a column
      * @param j
      * @return 
      */
     ArrayList<String> get_states_column(int j) {
         ArrayList<String> temp=new ArrayList<String>();
         for (int i=0; i<this.ntax;i++) {
             String d=this.char_matrix[i][j];
             for (int k=0; k<d.length();k++) {
                 char e=d.charAt(k);
                 if (!temp.contains(""+e)) temp.add(""+e);
             }
         }
         return temp;
     }

//     /**
//      * This handle the creation of bipartite graph
//      * @param node1
//      * @param node2 
//      */
//     public void compute_bipartition(node node1, node node2) {
//         BitVector ids=util.intersection_Bit(node1, node2);
//         int bipartite_id=bipartite_index;
//         String id=ids.toString();
//         if (!ids.equals(new BitVector(ntax))) {
//             if (bipartite_node_id.containsKey(ids)) {
//                   bipartite_id=bipartite_node_id.get(ids);
//             } else {
//                  bipartite_node_id.put(id,bipartite_index);
//                  bipartite_index++;
//             }
//         } 
//                    ArrayList<Integer> tmp=util.intersectBitResult(node1.partition,node2.partition);
//                     int total=tmp.size();
//                         //--Get the type here
//                        int type=4; //--default
//                        if (total==node1.total_taxa&&total==node2.total_taxa) {
//                            type=1;
//                        } else 
//                        if (total<node2.total_taxa&&total==node1.total_taxa) {
//                            type=2;
//                        } else 
//                        if (total<node1.total_taxa&&total==node2.total_taxa) {
//                            type=5;
//                        } else if (total>0) type=3;
//                       
//                        if (type!=4) {
//                            int source_index=node1.id;
//                            int dest_index=node2.id;
//                            if (type==5) {
//                                type=2;
//                                int tt=source_index;
//                                source_index=dest_index;
//                                dest_index=tt;                                
//                            } 
//                           
//                                if (bipartite_type==0||type==bipartite_type) {
//                                    output_biparition_complete.println(bipartite_id+"\t"+source_index);
//                                    output_biparition_complete.println(bipartite_id+"\t"+dest_index);
//                                } 
//                            
//                        } //--End not 4         
//     }
     
     
     public void compute_persistance_and_bipartition(node node1, node node2) {
                        BitVector ids=util.intersection_Bit(node1, node2);
                        int bipartite_id=bipartite_index;
                        String id=ids.toString();
                        if (!ids.equals(new BitVector(ntax))) {
                            if (bipartite_node_id.containsKey(id)) {
                                  bipartite_id=bipartite_node_id.get(id);
                            } else {
                                 bipartite_node_id.put(id,bipartite_index);
                                 bipartite_index++;
                            }
                        } 
                        
                        ArrayList<Integer> tmp=util.intersectBitResult(node1.partition,node2.partition);
                        int total=tmp.size();
                        //--Get the type here
                        int type=4; //--defaut
                        if (total==0) {
                            type=4;
                        } else if (total==node1.total_taxa&&total==node2.total_taxa) {
                            type=1;
                        } else 
                        if (total<node2.total_taxa&&total==node1.total_taxa) {
                            type=2;
                        } else 
                        if (total<node1.total_taxa&&total==node2.total_taxa) {
                            type=5;
                        } else 
                        if (total>0) {
                            type=3;
                        }
                       
                        if (type!=4) {
                            int source_index=node1.id;
                            int dest_index=node2.id;
                            if (type==5) {
                                type=2;
                                int tt=source_index;
                                source_index=dest_index;
                                dest_index=tt;                                
                            } 
                           if (bipartite) {
                            // Output bipartition
                                    output_biparition_complete.println(bipartite_id+"\t"+source_index+"\t"+type);
                                    output_biparition_complete.println(bipartite_id+"\t"+dest_index+"\t"+type);
                                    
                                    if (type==1) {
                                          output_biparition_1.println(bipartite_id+"\t"+source_index);
                                          output_biparition_1.println(bipartite_id+"\t"+dest_index);
                                    }
                                    if (type==2) {
                                          output_biparition_2.println(bipartite_id+"\t"+source_index);
                                          output_biparition_2.println(bipartite_id+"\t"+dest_index);
                                    }
                                    if (type==3) {
                                          output_biparition_3.println(bipartite_id+"\t"+source_index);
                                          output_biparition_3.println(bipartite_id+"\t"+dest_index);
                                    }
                           }
                            // Output persistance
                             int current_edge=total_edge;                             
                             total_edge++;
                             node1.edgecount++;
                             node2.edgecount++;
                             if (node1.id==source_index) node1.out_edgecount++;
                             if (node2.id==source_index) node2.out_edgecount++;
                             if (node1.id==dest_index) node1.in_edgecount++;
                             if (node2.id==dest_index) node2.in_edgecount++;
                            src_edge[current_edge]=source_index;
                            dest_edge[current_edge]=dest_index;
                            type_edge[current_edge]=type;
                            taxa_edge[current_edge]=total;
                           //count_edge[current_edge]=estimate_link_likelyhood(node1,node2);
                            nodes.set(node1.id, node1);
                            nodes.set(node2.id, node2);
                            
                                //--Statistics
                                total_type0++;
                                switch( type_edge[current_edge]) {
                                    case 1: total_type1++;break;
                                    case 2: total_type2++;break;
                                    case 3: total_type3++;break;
                                }
                                 node_id_type.get(type_edge[current_edge]).put(source_index, type);
                                 node_id_type.get(type_edge[current_edge]).put(dest_index, type);
                                 node_id_type.get(0).put(source_index, type);
                                 node_id_type.get(0).put(dest_index, type);                                 
        }

    }   
     
    
    /**
     * Return the taxa_id corresponding with the List
     * but with id starting at 1
     * @param ids_one
     * @return 
     */    
    String get_taxa(ArrayList<Integer> ids_one) {
        String str="";
        for (Integer id:ids_one) {
            str+=this.label.get(id-1)+",";
        }
        if (str.length()>2) str=str.substring(0,str.length()-1);
        return str;
    }
    
    String get_taxa(BitVector s) {
        String str="";
        for (int i=0; i<ntax;i++) {
            if (s.getBool(i)) str+=this.label.get(i)+",";
        }
        if (str.length()>1) str=str.substring(0,str.length()-1);
        return str;
    }
    
     String get_taxa(String s) {
        String str="";
        int index=0;
        for (int i=s.length()-1; i>-1;i--) {
            char c=s.charAt(i);
            if (c=='0') {
              index++;  
            } else if (c=='1'){
                 str+=this.label.get(index)+",";
                index++;
            } else {
                //Nothing
            }           
        }
        if (str.length()>1) str=str.substring(0,str.length()-1);
        return str;
    }
    
    BitVector get_taxa_bit(ArrayList<Integer> ids_one) {
        BitVector b=new BitVector(this.ntax);
        for (int i:ids_one) b.setBool(i-1, true);
        return b;
    }
    
    String get_taxa_id(ArrayList<Integer> ids_one) {
        String str="";
        for (Integer id:ids_one) {
            str+=id+",";
        }
        if (str.length()>1) str=str.substring(0,str.length()-1);
        return str;
    }
    
    //-- TO DO: only generate the node once
    public int estimate_link_likelyhood(node node1, node node2) {
        int total=0;
        int len=node1.partition.size();
        //1. Take the possible sate
        if (node1.multistate==1&&node2.multistate==1) return 0;
        // 2. compute each bitvector to estimate the number of time a link can be done.
        // ... we don't care now what kind of link...
        ArrayList<BitVector>node1_bit= new  ArrayList<BitVector>();
        ArrayList<BitVector>node2_bit= new  ArrayList<BitVector>();
        if (node1.multistate==1) {
            node1_bit.add(node1.partition);
        } else {             
            
              ArrayList<String> s1=util.combinations(get_column(node1.index-1));
              for (String p:s1) {
                BitVector tmp= new BitVector(len);
                 for (int j=0; j<len;j++) {
                   String p2=""+p.charAt(j);
                   if (p2.equals(node1.state_matrix)) {                      
                      tmp.setBool(j, true);
                  }
                }
                node1_bit.add(tmp);
              }
        }
        
         if (node2.multistate==1) {
        node2_bit.add(node2.partition);
        } else {
            ArrayList<String> s2=util.combinations(get_column(node2.index-1));
              for (String p:s2) {
                BitVector tmp= new BitVector(len);
                 for (int j=0; j<len;j++) {
                   String p2=""+p.charAt(j);
                   if (p2.equals(node2.state_matrix)) {                      
                      tmp.setBool(j, true);
                  }
                }
                node2_bit.add(tmp);
              }              
        }
       ///--Now, calculate the correspondance
       int total_possible=(node1_bit.size()*node2_bit.size());
//         for (BitVector b1:node1_bit) 
//           for (BitVector b2:node2_bit) 
//            if (util.intersectBitResult(b1,b2).size()>0) total++;
//            total=(100*total)/total_possible; //--compute the score
        return total;
    }
    
    /**
     * Given a node with a state s, evaluate the number of bitvector state
     * @param node1
     * @return 
     */
    public int get_node_multistate(node node1) {         
         int total=1;
          for (int i=0; i<this.ntax;i++) {
             String d=this.char_matrix[i][node1.index-1];
               if (d.length()>1) { 
                   if (d.indexOf(node1.state_matrix)>-1) {
                     total*=2; //--because its a binary state (either we have, or not)
                 }                    
               }            
         }         
         return total;
    }

    void load_charstate(String filename) {
        if (filename.isEmpty()) return;
        System.out.println("Loading external char.-states:"+filename);
        ArrayList<String> tmp=util.loadStrings(filename);
        if (tmp.size()==0) {
            System.out.println("...failed");
            return;
        }
       
        this.statelabels.clear();
        this.charlabels.clear();
        // For now, we expect the state labels in orer
        int last_ch=-1;
        String current_state="";
        ArrayList<String> tmp_ch=new ArrayList<String>();        
        for (String stri:tmp) {
            String[] s=stri.split("\t");
            if (s.length>=4) {
                int ch=Integer.valueOf(s[0]);
                //int state2=Integer.valueOf(s[1]);
                String states2=s[2].trim();
                String chs=s[3].trim();
                
                if (ch!=last_ch) {
                    if (last_ch!=-1) {
                        //System.out.println(last_ch+" "+current_state+tmp_ch);
                        this.charlabels.add(current_state);
                        this.statelabels.add(tmp_ch);
                    }
                    tmp_ch=new ArrayList<String>();
                    last_ch=ch;
                    current_state=states2;
                }
                tmp_ch.add(chs);
            }  else {
                System.out.println("Error (not enough tabulation): "+stri);
            }  
        }
        this.charlabels.add(current_state);
        this.statelabels.add(tmp_ch);       
//        for (int i=0; i<charlabels.size();i++) {
//            System.out.println((i+1)+" "+charlabels.get(i)+" "+statelabels.get(i));
//        }
    }
    
    
    
    public void analyse(String filename, int type) {
        graph g=new graph();
        g.load_graph(filename, false);
        //--Print out statistics 
        String filename_out=filename+".summary.txt";
        
    }
    
} //End dataset
