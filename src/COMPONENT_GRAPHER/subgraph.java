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

/**
 * This create a subgraph of the main graph so that we can compute local attributes
 * @author Etienne Lord, Jananan Pathmanathan
 * @since October/November 2015
 */
public class subgraph extends graph {
    
    HashMap<Integer,Integer> oldid_to_newid= new HashMap<Integer,Integer>();
    /**
    * This construct a subgraph centered around node s
    */
    public subgraph(graph g,int s) {
        this.total_nodes=g.get_adj(s).size();
        for (int n:g.get_adj(s).keySet()) {
            // Reorder node
            Integer id=oldid_to_newid.get(n);
            if (id==null) {
                id=oldid_to_newid.size();
                oldid_to_newid.put(n,id);
            } 
        }
          for (int n:g.get_adj(s).keySet()) {
            Integer id=oldid_to_newid.get(n);
            HashMap<Integer,Boolean> tmp=new HashMap<Integer,Boolean>();
            for (int w:g.get_adj(n).keySet()) if (oldid_to_newid.get(w)!=null&&w!=s) tmp.put(oldid_to_newid.get(w), Boolean.TRUE);
            this.adjlist.put(id, tmp);
        }
    }
    
}

