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

/**
 *.Class representing a character-state from the input matrix
 * @author Etienne Lord, Jananan Pathmanathan
 * @since October/November 2015
 */
public class state {
    
    ////////////////////////////////////////////////////////////////////////////
    /// VARIABLES
    int state_id=0;  //--ID for this state
    int pos_i=0;     //--Position in matrix
    int pos_j=0;     //--Position in matrix
    ArrayList<String> states=new ArrayList<String>(); //--State as en ArrayList
    String state=""; //--State at this position (e.g. ABC)    
}
