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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import umontreal.iro.lecuyer.util.BitVector;
import java.io.*;
import java.util.Calendar;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Loading, intersection detection, etc.
 * @author Etienne Lord, Jananan Pathmanathan
 * @since October/November 2015
 */
public class util {
    
   public static HashMap<Integer,Long> bionomial_cache=new  HashMap<Integer,Long>();
   public static Pattern isNumber=Pattern.compile("[0-9]{1,}(.[0-9]*])*");
   public PrintWriter pw; //--global object for fast output
   
   public static ArrayList<Integer> union(ArrayList<Integer> list1, ArrayList<Integer> list2) {
      ArrayList<Integer> list = new ArrayList<Integer>();
       for (Integer t : list1) {
            if(!list.contains(t)) {
                list.add(t);
            }
        }
        for (Integer t : list2) {
            if(!list.contains(t)) {
                list.add(t);
            }
        }
       return list;
   }
        
     public static <T> List<T> intersection(List<T> list1, List<T> list2) {
        List<T> list = new ArrayList<T>();

        for (T t : list1) {
            if(list2.contains(t)) {
                list.add(t);
            }
        }

        return list;
    }
     
     /**
      * This will determine the intersection type of two list
      * Note: undirected 1 vs 2
      * @param list1
      * @param list2
      * Type 1: total 1,2,3 vs 1,2,3
      * Type 2: partial 1,2,3 vs 1,2,3,4,5 ->
      * Type 5: partial  1,2,3,4,5 vs 1,2,3 <-
      * Type 3: partial 1,2,3 vs 2,3,4
      * Type 4: disjoint 1,2,3 vs 4,5,6
      */
     public static <T> int intersection_type(List<T> list1, List<T> list2) {
         List<T> tmp=intersection(list1,list2);
         if (tmp.size()==list2.size()&&tmp.size()==list1.size()) return 1;
         if (tmp.size()<list2.size()&&tmp.size()==list1.size()) return 2;
         if (tmp.size()<list1.size()&&tmp.size()==list2.size()) return 5;
         if (tmp.size()>0) return 3;
         return 4; 
     }
    
     public static ArrayList<ArrayList<Integer>> get_partition(String data) { 	
	// Note, the first element is element 1 in the returned
	// data.
	
	ArrayList<ArrayList<Integer> > combinations=new ArrayList<ArrayList<Integer>>();
	HashMap<Character, ArrayList<Integer> > d=new HashMap<Character, ArrayList<Integer> >();	
	
	for (int j=0; j<(int)data.length();j++) {
		//Note, we skip not defined character for now
		if (data.charAt(j)!='?'&&data.charAt(j)!='-'&&data.charAt(j)!='*') {
			ArrayList<Integer> tmp=new ArrayList<Integer>();
                        if (d.containsKey(data.charAt(j))) {
                            tmp= d.get(data.charAt(j));				
			}
                        tmp.add(j+1);
                        d.put(data.charAt(j),tmp);
		}
	}
       
         HashMap<Integer,Character> keys=new HashMap<Integer,Character>();
        for (Character k:d.keySet()) {
            Integer v=d.get(k).get(0);
            keys.put(v,k);
        }
        
       for (Integer i:keys.keySet()) {
           combinations.add(d.get(keys.get(i)));
       }
	return combinations;
     }
  
      public static ArrayList<ArrayList<Integer>> get_partition(ArrayList<String> data) { 	
	// Note, the first element is element 1 in the returned
	// data.	
	ArrayList<ArrayList<Integer> > combinations=new ArrayList<ArrayList<Integer> > ();
	HashMap<String, ArrayList<Integer> > d=new HashMap<String, ArrayList<Integer> >();	
	
	for (int j=0; j<(int)data.size();j++) {
		//Note, we skip not defined character for now
		if (!data.get(j).equals("?")&&!data.get(j).equals("-")&&!data.get(j).equals("*")) {
			ArrayList<Integer> tmp=new ArrayList<Integer>();
                        if (d.containsKey(data.get(j))) {
                            tmp= d.get(data.get(j));				
			}
                        tmp.add(j+1);
                        d.put(data.get(j),tmp);
		}
	}
       
         HashMap<Integer,String> keys=new  HashMap<Integer,String> ();
        for (String k:d.keySet()) {
            Integer v=d.get(k).get(0);
            keys.put(v,k);
        }
        
       for (Integer i:keys.keySet()) {
           combinations.add(d.get(keys.get(i)));
       }
	return combinations;
     }
      
      /**
     * This create from a multiple data of String e.g.
     * 12,3,45,6 
     * The multiple combinations:
     * 1346
     * 2346
     * 1356
     * 2356
     * @param data
     * @return All the encoded string
     */
    public static  ArrayList<String> combinations(ArrayList<String> data) { 	
	 ArrayList<String> tmp=new  ArrayList<String>();        
	tmp.add(""); 
	for(int i=0;i<data.size();i++) {
		String car=data.get(i);
		
                if (car.length()>1) {
				ArrayList<String> v=new ArrayList<String>();
                                for (int j=0; j<car.length();j++) v.add(""+car.charAt(j));				
				// Duplicate vector here
				int pos=v.size()-1;
				int pos2=tmp.size();
				for (int j=0; j<pos;j++) {
					for (int k=0; k<pos2;k++) {
					    String ctmp=new String(tmp.get(k));
                                            tmp.add(ctmp);	
					}
				}	
				int pos_i=0;			
				for (int j=0; j<v.size();j++) {
					for (int k=0; k<pos2;k++) {
                                               String s=tmp.get(pos_i);
						tmp.remove(pos_i);
                                               s+=v.get(j);
                                                tmp.add(pos_i++, s);
					}
				}                                		
		} else {
                       int l=tmp.size();
			for (int j=0; j<l;j++) {
			    String s=tmp.get(j);
                            tmp.remove(j);
                            s+=car;    
                            tmp.add(j,s);
			}				
		}		
	}       
	return tmp;
 }
     /**
      *  Return the binomial of number n by k
      * @param n
      * @param k
      * @return 
      */
     private static long binomial(int n, int k)
    {
       if (bionomial_cache.containsKey(n)) {
           return bionomial_cache.get(n);
       }
        if (k>n-k)
            k=n-k;
 
        long b=1;
        for (int i=1, m=n; i<=k; i++, m--)
            b=b*m/i;
        bionomial_cache.put(n, b);
        return b;
    }
     
    /**
     * Rand Index of two BitVector
     * @param b1
     * @param b2
     * @return the RandIndex value
     */
    static double RandIndexBit(BitVector b1, BitVector b2)  {
        double a=0;
        double b=0;
        for (int i=0; i<b1.size();i++) {
            boolean t1=b1.getBool(i);
            boolean t2=b2.getBool(i);
            for (int j=0;j<b1.size();j++) {
                if (i<j) {
                    boolean t3=b1.getBool(j);
                    boolean t4=b2.getBool(j);
                    if (t1==t3&&t2==t4) a++;
                    if (t1!=t3&&t2!=t4) b++;
                }
            }
        }
        return (a+b)/binomial(b1.size(),2);
    }
    
    /**
     * Return the Number of common bits between bit vector
     * @param b1
     * @param b2
     * @return the number of common bits
     */
    public static int intersectBit(BitVector b1, BitVector b2) {
        int total=0;
        BitVector b3=b1.and(b2);
        for (int i=0;i<b3.size();i++) 
            if (b3.getBool(i)) total++;
        return total;
    }
    
   
    /**
     * Return the intersecting values of two BitVector
     * @param b1
     * @param b2
     * @return the position with index starting at 1
     */
     public static ArrayList<Integer> intersectBitResult(BitVector b1, BitVector b2) {
        ArrayList<Integer> tmp=new ArrayList<Integer>();        
        BitVector b3=b1.and(b2);
        for (int i=0;i<b3.size();i++) 
            if (b3.getBool(i)) tmp.add(i+1);
        return tmp;
    }
    
     /**
      * Return the number of set (true) bit in a BitVector
      * @param b1
      * @return the value
      */
    public static int total_bitset(BitVector b1) {
        int total=0;
         for (int i=0;i<b1.size();i++) 
            if (b1.getBool(i)) total++;        
        return total;
    }
   
    /**
     * Return a new BitVector contained in nodes from the AND operation
     * @param n1 
     * @param n2
     * @return 
     */
      public static BitVector intersection_Bit(node n1, node n2) {
           BitVector b3=n1.partition.and(n2.partition);
          return b3;
      }
      
      /**
       * Open a PrintWriter for easy file output
       * @param filename 
       */
     public void open(String filename) {
         try {
             pw=new PrintWriter(new FileWriter(new File(filename)));
         } catch(Exception e) {}
     } 
     
     /**
      * Print stri to the open PrintWriter
      * @param stri 
      */
     public void println(String stri) {
         try {
             pw.println(stri);
         } catch(Exception e) {}
     }
    
     /**
      * Close the PrintWriter
      */
     public void close() {
         try {
             pw.flush();
             pw.close();
         } catch(Exception e) {}
     }
     
     /**
      * Load a filename into an ArrayList of Strings
      * @param filename
      * @return an ArrayList<String>
      */
     public static ArrayList<String> loadStrings(String filename) {
         ArrayList<String> tmp=new ArrayList<String>();
         try {
             BufferedReader br=new BufferedReader(new FileReader(new File(filename)));
             while(br.ready()){
                 tmp.add(br.readLine());
             }
             br.close();
          } catch(Exception e) {}
         
         return tmp;
     }
     
      /**
     * Test if a directory exists
     * @param filename
     * @return true if file Exists
     */
    public static boolean DirExists(String filename) {
             File f = new File(filename);
             if (f==null||!f.isDirectory()) return false;
             return f.exists();
    }
    
    /**
     * Create the directory or the directory arborescence 
     * @param directory
     * @return true if succes, false for any other reasons.. (dir exists...)
     */
    public static boolean CreateDir(String directory) {
        try {
            File f=new File(directory);            
            if (DirExists(directory)) return false;
            return (f.mkdirs());           
        } catch(Exception e) {
            return false;
        }
    }
    
    /**
     * Using Regex, find if we have a number in the String
     * @param data
     * @return 
     */
    public static boolean isNumber(String data) {
        Matcher m=isNumber.matcher(data);
        return m.find();
    }
    
    /**
     * Convert an array of String into one String
     * @param data
     * @return 
     */
    public static String string(String[] data) {
        String st="";
        for (String s:data) st+=s+" ";     
        return st;
    }
        
     /**
     * This return a time code JJMMYYMMSS used by the hashcode
     * @return
     */
    public static int returnTimeCode() {
        Calendar today=Calendar.getInstance();
        int[] tt={0,0,0,0,0,0};
        tt[0]=today.get(Calendar.DAY_OF_MONTH);
        tt[1]=today.get(Calendar.MONTH);
        tt[2]=today.get(Calendar.HOUR);
        tt[3]=today.get(Calendar.MINUTE);
        return (tt[0]*1000000+tt[1]*10000+tt[2]*100+tt[3]);
    }
    
    /**
        * Transform a number in ms to its d,h,m,s representation
        * @param ms (system millisecond)
        * @return a String representation
        */
       public static String msToString(long ms) {
            long d, h, m, s;
            d=ms / 86400000;
            ms=ms-(d*86400000);
            h = ms / 3600000;
            ms=ms-(h*3600000);
            m=ms/60000;
            ms=ms-(m*60000);
            s=ms/1000;
            ms=ms-(s*1000);
            String str="";
            if (d>2000) return ms+" ms "; //hack
            if (d!=0) return d+"d "+h+"h "+m+"m "+s+"s ";
            if (h!=0) return h+"h "+m+"m "+s+"s ";
            if (m!=0) return m+"m "+s+"s ";
            if (s!=0) return s+"s "+ms+"ms ";
            return ms+"ms ";
       }
       
    /**
     * List the files in directory and return filename with the full path
     * @return
     */
    public static ArrayList<String> listDirWithFullPath(String dir) {
        ArrayList<String> tmp=new ArrayList<String>();
        try {
            File f=new File(dir);
            if (!f.isDirectory()) {
                dir+=File.separator;
                f=new File(dir);
            }
            for (String filename:f.list()) tmp.add(f.getAbsolutePath()+File.separator+filename);
        } catch(Exception e) {}
        return tmp;
    }
    
       
       
       /**
        * Recursive function to list a directory 
        * Note: Use the listDirFullPath function
        * @param path
        * @return an ArrayList of the full file path
        */
       public static ArrayList<String> recursive_list_dir(String path) {           
       ArrayList<String>filenames = new ArrayList<String>();
       ArrayList<String>  filename_tmp=listDirWithFullPath(path);        
        for (String fi:filename_tmp)        
        {
            File f=new File(fi);             
            if (f.isDirectory()) {
                filenames.addAll(recursive_list_dir(f.getAbsolutePath()));                
            } else {
                filenames.add(fi);
            }                     
        }
        return filenames;
    }

    /**
     * This filter a list of files according to a regex
     * @param filenames
     * @param regex
     * @return 
     */   
     public static ArrayList<String>filter(ArrayList<String> filenames, String regex) {
         ArrayList<String>tmp=new ArrayList<String>();  
         Pattern p=Pattern.compile(regex);    
         for (String s:filenames) {
           Matcher m=p.matcher(s); 
            if (m.find()) tmp.add(s); 
         }
         return tmp;
     }
   
}
