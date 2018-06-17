package grouprepresentation;

import java.io.FileWriter;

/** Functions for the re-implementation of Faulon's paper have been developed.
 * 
 * 1. The main idea is how to check whether graph is ideal (R) or not.
 * 
 * 2. If graph is not ideal, needs be made an ideal graph.
 * 
 * 3. In (R), permutations were used for equivalence classes. Here, Canon.symmetry
 * class of CDK is implemented. The class calculates the symmetric equivalences of vertices
 * based on the permutations.
 * 
 * 4. For equivalence class classification and generation, sub functions are coded.
 * 
 * 5. The initial graph is extended from the equivalence classes' representatives.
 * In other words, no need to extend all the vertices but just their unique equivalence
 * in the graph.
 * 
 * For detailed explanation, please take a look at the referenced paper.
 * 
 * Reference (R): Faulon, Jean Loup. "On using graph-equivalent classes for the structure
 * elucidation of large molecules." Journal of chemical information and computer sciences 
 * 32.4 (1992): 338-348.
 */


//Functions might be removed: bondcounter, molopens,
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.SpanningTree;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.graph.invariant.EquivalentClassPartitioner;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.BondManipulator;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;

import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor;

public class Functions  {
	static SDFWriter outFile;
	private static Functions instance;
	private static SaturationChecker saturation;
    
	public ArrayList<IAtomContainer> aclist;
 
	private Functions() throws IOException {
	    aclist = new ArrayList<IAtomContainer>();
	    outFile = new SDFWriter(new FileWriter("C:\\Users\\MAY\\Desktop\\sonuc.sdf"));
	}
		
	public static Functions getInstance() throws IOException {
	    if (instance == null) instance = new Functions();
	    return instance;
	  }

	
/**
 * These are the general functions used for the basics.
 * 
 */

	public static ListMultimap<Long, Long> repclean(ListMultimap<Long, Long> multmap){
		
		//To clean the reputation like 1={5} 5={1} for symmetry classes
		for(long key: multmap.keys()){
			for( long l:multmap.get(key)){ //The classes of the values are checked.
					if(key!=l && multmap.get(l).contains(key) &&multmap.get(l).size()!=1){ // key!=l is for 3={3}, removing the same one from one of the classes. no more 1={5} 5={1} for symmetry classes
						multmap.get(l).remove(key);
				}
			}
		}
		return multmap;
	}
	
	//To get the index of long array
	public static int indexOf(long[] array, long target) {
		for (int i = 0; i < array.length; i++) {
	        if (target==array[i]){
	            return i;
	        }
	    }
	    return -1;
	}

	//Store all the atoms' connected bonds in an array.
	public static int[] bondcounter (IAtomContainer mol) throws CloneNotSupportedException, CDKException {
		int [] bondnumbers= new int [mol.getAtomCount()];
		
		//Count the interactions for all the vertices and store in the array
		for(IAtom atom: mol.atoms()){
			bondnumbers[mol.indexOf(atom)]=mol.getConnectedBondsList(atom).size();
		}
				
		return bondnumbers;
	}
	
	// Molecule depiction
	public static void depict(IAtomContainer mol, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depict = new DepictionGenerator();
		depict.withSize(200, 250).withAtomColors().depict(mol).writeTo(path); 
		
	}
			
	//Counting hydrogens in the atomcontainer
	public static List<IAtom> hydcount(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		List<IAtom> hydr = new ArrayList<IAtom>();
		for(IAtom atom: mol.atoms()){
			if(atom.getSymbol()=="H"){
				hydr.add(atom);
			}
		}
	return hydr;
	}
	
			
	// Removing the hydrogens from the atomcontainer
	public static void removehyd( IAtomContainer mol, List<IAtom> hydr) throws CloneNotSupportedException, CDKException, IOException {
		for(IAtom atom: hydr){
			mol.removeAtom(atom);
		}
	}
			
	//Counting open bond sites
	public static int molopens(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		int open=0;
		for(IAtom atom: mol.atoms()){
			if(atom.getSymbol()=="H"){
				open--;
			}
			else{
				open += valences.get(atom.getSymbol());
			}
		}
		return open;
	}
		
	// Distinguishing the cyclic and non-cyclic atoms in a molecule
	public static List<IAtom> ringcheck (IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException {
		SpanningTree   tree   = new SpanningTree(mol);
		IAtomContainer cyclic = tree.getCyclicFragmentsContainer();
		List<IAtom> inring = new ArrayList<IAtom>();
			
		// Detect the atoms appeared in ring structures.
		for(IAtom atom : mol.atoms()){
			if(cyclic.contains(atom)){
			    inring.add(atom);
			}
		}
		return inring;
	}
		
	//This functions classifies the atoms with their sym class representatives; then returns the multmap.
	public static ListMultimap<String, Long> atom2symc (IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException {
		long[] sym= Functions.canonsym(mol);
	    ListMultimap<String, Long> multmap = ArrayListMultimap.create();
	    for(int s=0; s<mol.getAtomCount();s++){ //The atoms are classified with their symbols and equivalent sym class.
	    	if(!multmap.get(mol.getAtom(s).getSymbol()).contains(sym[s])){//For extending the atoms, type and sym equivalences pairs are listed.
	    		multmap.put(mol.getAtom(s).getSymbol(), sym[s]);
	        }
	    }
	    return multmap;
	}
	
	//Rather then returning the sym map values of the element types, return their possible extension indices based on symmetry values.
	public static ListMultimap<String, Integer> atom2symind (IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException {
		long[] sym= Functions.canonsym(mol);
		ListMultimap<String, Long> atom2sym = atom2symc(mol);
	    ListMultimap<String, Integer> atom2symindex = ArrayListMultimap.create();

	    for(String key: atom2sym.keys()){ 
	    	for(Long val:atom2sym.get(key)){
	    		if(!atom2symindex.get(key).contains(Functions.indexOf(sym, val))){//For extending the atoms, type and sym equivalences pairs are listed.
	    			atom2symindex.put(key, Functions.indexOf(sym, val));
		        }
	    	}
	    }
	    return atom2symindex;
	}
	
	//This function is the original version of equivalence class classification considering also the implicit hydrogens
	public static ListMultimap<String, Integer> atom2impclass(IAtomContainer acontainer) throws CloneNotSupportedException, CDKException, IOException {
		ListMultimap<String, Integer> classs = ArrayListMultimap.create();
		//Long.valueOf(sym[i]).intValue();
		for(int i=0; i<acontainer.getAtomCount();i++){
			if(Functions.satcheck(acontainer, i)==true){
				classs.put(acontainer.getAtom(i).getSymbol()+acontainer.getAtom(i).getImplicitHydrogenCount(), i);
			}
		}
		
		return classs;
	}
	
	//This function is the original version of equivalence class classification considering also the implicit hydrogens
	public static ListMultimap<String, Integer> atom2impclass2(IAtomContainer acontainer) throws CloneNotSupportedException, CDKException, IOException {
		ListMultimap<String, Integer> classs = ArrayListMultimap.create();
		long[] sym=Functions.canonsym(acontainer);
		//System.out.println("sym"+ " "+Arrays.toString(sym));
		//Long.valueOf(sym[i]).intValue();
		for(int i=0; i<acontainer.getAtomCount();i++){
			if(Functions.satcheck(acontainer, i)==true){
			//System.out.println("tom2imp"+ " "+i+ " "+acontainer.getAtom(i).getSymbol()+" "+acontainer.getAtom(i).getImplicitHydrogenCount()+" "+Long.valueOf(sym[i]).intValue());
				classs.put(acontainer.getAtom(i).getSymbol()+acontainer.getAtom(i).getImplicitHydrogenCount()+Long.valueOf(sym[i]).intValue(), i);
			}
		}
			
		return classs;
	}
	//ListmultiMap includes the key duplicates, so this is just for cleaning the duplicates.
	public static ArrayList<String> cleanatom2impclass(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		ListMultimap<String, Integer> map = Functions.atom2impclass(mol);
		ArrayList<String> keys= new ArrayList<String>();
		for(String l:map.keys()){
			if(keys.contains((String)l)) continue;
		    	keys.add(l);
		    }
		return keys;
	}
	
	//ListmultiMap includes the key duplicates, so this is just for cleaning the duplicates.
	public static ArrayList<String> cleanatom2impclass2(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		ListMultimap<String, Integer> map = Functions.atom2impclass2(mol);
		ArrayList<String> keys= new ArrayList<String>();
		for(String l:map.keys()){
			if(keys.contains((String)l)) continue;
				keys.add(l);
			}
			return keys;
	}
	//This function is the original version of equivalence class classification considering also the implicit hydrogens
	public static ListMultimap<Integer, Integer> atom2impclasshashcode(IAtomContainer acontainer) {
		long[] sym= Functions.canonsym(acontainer);
		ListMultimap<Integer, Integer> classs = ArrayListMultimap.create();
		for(int i=0; i<acontainer.getAtomCount();i++){
			classs.put(acontainer.getAtom(i).hashCode(), Long.valueOf(sym[i]).intValue());
		}
		return classs;
	}
	
	//This function is the original version of equivalence class classification considering also the implicit hydrogens
	public static List<Object> classrep(IAtomContainer acontainer) throws CloneNotSupportedException, CDKException, IOException {
		List<Object> list = new ArrayList<Object>();
    	ListMultimap<String, Integer> map= Functions.atom2impclass(acontainer);
    	for(String key:Functions.cleanatom2impclass(acontainer)){
    		list.add(key);
    		Collections.sort(map.get(key));
    		for(Integer val: map.get(key)){
    			list.add(val);
    		}
    	}
		return list;
	}
	
	//This function is the original version of equivalence class classification considering also the implicit hydrogens
	public static List<Object> classrep2(IAtomContainer acontainer) throws CloneNotSupportedException, CDKException, IOException {
		List<Object> list = new ArrayList<Object>();
	    ListMultimap<String, Integer> map= Functions.atom2impclass2(acontainer);
	    //System.out.println("map"+ " "+map);
	    for(String key:Functions.cleanatom2impclass2(acontainer)){
	    	list.add(key);
	    	Collections.sort(map.get(key));
	    	for(Integer val: map.get(key)){
	    		list.add(val);
	    	}
	    }
		return list;
	}
	
	//This function is the identical graph checker comparing all the equivalence classes whether they are identical or not. If identical, then the graphs are isomorphic.
	public static boolean identical_class (IAtomContainer ac1, IAtomContainer ac2) throws CloneNotSupportedException, CDKException, IOException{
		boolean check =true;
        if(Functions.classrep2(ac1).size()!=Functions.classrep2(ac2).size()){
        	check=false;
        }else{
        	for(int i=0;i<Functions.classrep2(ac1).size();i++){
            	if(!Functions.classrep2(ac1).get(i).equals(Functions.classrep2(ac2).get(i))){
            		check=false;
            		break;
            	}
            }
        }
		return check;
	}
	
/**
* These functions are for cleaning the interactions
* of the same atom types, if needed.
*/

	//For collecting the atoms of the same type.
	//stype means same type.
	public static List<IAtom> sametypes(IAtomContainer mol, String sym){
		List<IAtom> same = new ArrayList<IAtom>();
		for(IAtom atom: mol.atoms()){
			if(atom.getSymbol()==sym){
				same.add(atom);
			}
		}
		return same;
	}
	
	//Collects the same typed atoms' connected bonds.
	public static List<IBond> sbonds(IAtomContainer mol, String sym){
		List<IBond> sbonds = new ArrayList<IBond>();
		for(IAtom atom: sametypes(mol, sym)){
			for(IBond bond: mol.getConnectedBondsList(atom)){
				sbonds.add(bond);
			}
		}
		return sbonds;
	}
	
	//To make a molecule ideal, this fucntions removes the same atoms' bonds distrubing the ideality case.
	public static void removebonds(IAtomContainer mol, String sym){
		List<IBond> bonds=sbonds(mol,sym);
		for(IBond bond: bonds){
			mol.removeBond(bond);
		}
	}
	
	//No need to call all the previous functions separately like in this one but just left as an example.
	public static void sintremove( IAtomContainer mol, String sym){
		//List<IAtom> stype= Functions.sametypes(mol, sym);
		//List<IBond> sbonds= Functions.sbonds(mol, sym);
		Functions.removebonds(mol, sym);
	}
	
/**
* These functions are for symmetry based equivalence class classification.
* 
*/
	
	//Calculates the canon symmetry array.
	public static long[] canonsym(IAtomContainer mol){
		int[][]	g = GraphUtil.toAdjList(mol);
    	long[] sym= Canon.symmetry(mol, g);
    	return sym;
	}
	
	//ListmultiMap includes the key duplicates, so this is just for cleaning the duplicates.
	public static ArrayList<Long> cleansymlistkey(IAtomContainer mol){
		ListMultimap<Long, Long> symlist = Functions.symlist(mol);
        ArrayList<Long> keys= new ArrayList<Long>();
        for(Long l:symlist.keys()){
        	if(keys.contains((long)l)) continue;
        	keys.add(l);
        }
        return keys;
	}

	
	// The function stores the symmetry classes with list of atoms.
	public static ListMultimap<Long, Long> symlist( IAtomContainer ac){
		long[] labels = Functions.canonsym(ac);
		ListMultimap<Long, Long> multmap = ArrayListMultimap.create();
		for(int i=0; i<labels.length;i++){ //The symmetry values of the atoms are stored based on their initial enumaration.
			multmap.get(labels[i]).add((long)i+1); //To start from 1.
			//multmap.put(labels[i], (long)i+1); //due to int long difference, it cannot realize the duplicates.
		}
		return multmap;
		}
	
/**
* For isomorphism check (CNI) interactions among the eclass elements
* and interaction with target are check with the following functions. 
*/
	
	//Checks inner interactions among class elements.
	public static boolean intinclass(IAtomContainer mol, Long key){
		ListMultimap<Long, Long> symlist=Functions.symlist(mol);
		Boolean check=true;
		for(int i=0;i<symlist.get(key).size();i++){
			for(int j=i+1;j<symlist.get(key).size();j++){
				if((mol.getBond(mol.getAtom(symlist.get(key).get(i).intValue()-1),mol.getAtom(symlist.get(key).get(j).intValue()-1)))!=null){
					check=false;
					break;
				}
			}
		}
		return check;
	}
	
	//Checks whether the class elements are interacted with the target or not.
	public static boolean intwithtarget(IAtomContainer mol, Long key) throws CloneNotSupportedException, CDKException, IOException{
		ListMultimap<Long, Long> symlist=Functions.symlist(mol);
		Boolean check=true;
		for(Long val:symlist.get(key)){
			if((mol.getBond(mol.getAtom(val.intValue()-1),mol.getAtom(Functions.classnext(mol, val.intValue()-1))))!=null){
				check=false;
				break;
			}
		}
		
		return check;
	}
	
	// The function stores the symmetry classes
	public static ListMultimap<Long, Long> symclass( IAtomContainer ac){
		long[] labels = Functions.canonsym(ac);
		ListMultimap<Long, Long> multmap = ArrayListMultimap.create();
		for(int i=0; i<labels.length;i++){ //The symmetry values of the atoms are stored based on their initial enumaration.
			multmap.put(labels[i], (long)i+1); //due to int long difference, it cannot realize the duplicates.
		}
		return repclean(multmap); //Cleaning the key reputations.
	}
		
	//Detects only the indices of the fixed non-symmetric atoms in sym array
	public static List<Integer> nonsymdetect(long[] sym){
		Map<Long, Long> freq = Arrays.stream(sym).boxed().collect(Collectors.groupingBy(Long::longValue, Collectors.counting()));
		List<Integer> indices = new ArrayList<Integer>(); //Indices of the non-sym atoms.
        for(long key:freq.keySet()){
        	if(freq.get(key)==1){ //key is sym value appeared once in the array.
        		indices.add(Functions.indexOf(sym, key));
        	}
        }
        return indices;
	}
	
/**
* These fucntions for checking whether the initial graph is ideal or not.
* If not, we can make it ideal. 
*/
	//TODO: General cases should be considered.
	//Takes molecule and its symmetry array from Canon.symmetry, and remove non symmetric ones' interactions not the atoms.
	public static void makeIdeal(IAtomContainer mol) throws NumberFormatException, CloneNotSupportedException, CDKException, IOException{
		long[] sym= Functions.canonsym(mol);
		List<Integer> indk= nonsymdetect(sym); //Nonsym indices are collected. Clean these types interactions. These are symmetry 
		for(int i:indk){ //For them, remove all the same atomtype interactions such as all the Cs.
			if(mol.getAtom(i).getSymbol()=="C" && Functions.ctypchecker(mol)==true){//Having single type of C means all the C atoms are connected to each other somehow.
				if(neigsym(mol,neigC(mol,mol.getAtom(i))).size()==0) continue; //If the other carbons also does not have any non-carbon interactions.
				else{
					for(String sm:neigsym(mol,neigC(mol,mol.getAtom(i)))){ //Else remove the non-carbon interactions of the mol, causing the sym destruction
						Functions.sintremove(mol, sm);
					}
					
				}
			}else{
				Functions.removeCbonds(mol);
				Functions.sintremove(mol, mol.getAtom(i).getSymbol());
			}
		}
	}
	
	//This function checks whether the carbons are the same type or not. If not, It removes all the carbon interactions for ideality case.
	//Functions.ctypchecker(mol)
	public static List<IBond> cbonds(IAtomContainer mol) throws NumberFormatException, CloneNotSupportedException, CDKException, IOException{
		List<IBond> cbnds= new ArrayList<IBond>();
		List<IAtom> cs= Functions.clist(mol);
		for(IAtom atm:cs){
			for(IAtom atm2:cs){
				if(mol.getBond(atm, atm2)!=null && !cbnds.contains(mol.getBond(atm, atm2))){
					cbnds.add(mol.getBond(atm, atm2));
				}
			}
        }
		return cbnds;
	}
	
	//Collecting the carbon atoms of the atom container to use in carbon bond collector. 
	public static List<IAtom> clist(IAtomContainer mol) throws NumberFormatException, CloneNotSupportedException, CDKException, IOException{
		List<IAtom> clst= new ArrayList<IAtom>();
		for(IAtom atm:mol.atoms()){
			if(atm.getSymbol()=="C"){
				clst.add(atm);
			}
        }
		return clst;
	}

	//Removing list of bonds.
	public static void removeCbonds(IAtomContainer mol) throws NumberFormatException, CloneNotSupportedException, CDKException, IOException{
		List<IBond> clst= Functions.cbonds(mol);
		for(IBond bond:clst){
			mol.removeBond(bond);
        }
	}
			
	//Collects the non carbon neighbors of the atom, if there are any.
	public static List<String> neigsym(IAtomContainer mol, IAtom atom) throws NumberFormatException, CloneNotSupportedException, CDKException, IOException{
		List<String> neig= new ArrayList<String>();
		for(IAtom atm:mol.getConnectedAtomsList(atom)){
        	if(!neig.contains(atm.getSymbol()) && atm.getSymbol()!="C"){
        		neig.add(atm.getSymbol());
        	}
        }
		return neig;
	}
	//Finds the neighbor carbon atom.
	public static IAtom neigC(IAtomContainer mol, IAtom atom) throws NumberFormatException, CloneNotSupportedException, CDKException, IOException{
		IAtom neigC= null;
		for(IAtom atm:mol.getConnectedAtomsList(atom)){
        	if(atm.getSymbol()=="C"){
        		neigC=atm;
        		break;
        	}
        }
		return neigC;
	}
	
	//This function basically check if all the atoms, based on symbols, are in the single sym class, then the bonds are also isomorphic.
	public static boolean idealcheck(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		boolean check= true;
		ListMultimap<String, Long> map=Functions.atom2symc(mol); // The sym representative based ranking.
		for(String key: map.keys()){ //Checking for symmetry destruction.
			if (map.get(key).size()!=1){ //Keys are the element symbols.
		        check=false;
		    }
		}
		return check;
	}
	
	//Check whether all the carbons are the same type in the molecule.
	public static boolean ctypchecker(IAtomContainer mol) throws NumberFormatException, CloneNotSupportedException, CDKException, IOException{
        CarbonTypesDescriptor desc= new CarbonTypesDescriptor(); 
        List<String> types = Arrays.asList(desc.calculate(mol).getValue().toString().split(","));
        boolean check= false;
        for(String key: types){
        	if(Functions.atomtypecount(mol, "C")==Integer.parseInt(key)){
        		check=true;
        	}
        } 
        return check;
	}
	
	//AtomCounter with CDK descriptors.
	public static int atomtypecount(IAtomContainer mol, String sym) throws CloneNotSupportedException, CDKException, IOException{
	// :) only create new
		AtomCountDescriptor descriptor =  new AtomCountDescriptor();
		Object[] params = {sym};
		descriptor.setParameters(params);   
		DescriptorValue value = descriptor.calculate(mol);    
		return (Integer.parseInt(value.getValue().toString())); //For comparisons, convert to int.
			
	}
			
/**
* These functions are specifically for extension process.
* 
*/
	//Get the total order of the connected bonds to an atom
	public static int ordsum(IAtomContainer mol, int i){
		int count=0;
		for(IBond bond: mol.getConnectedBondsList(mol.getAtom(i))){
        	count=count+bond.getOrder().numeric();
        }
		return count;
	}
	//Saturation check by comparing the maximum number of atoms can be attached and the currently attached one.
	public static boolean satcheck(IAtomContainer mol, int i)throws CloneNotSupportedException, CDKException, IOException{
		if ((mol.getAtom(i).getImplicitHydrogenCount()+ordsum(mol,i))>= (int)valences.get(mol.getAtom(i).getSymbol())){ //I tried with mol.getConnectedBondsCount(mol.getAtom(i))
			//System.out.println("satcheckvalues");
			//System.out.println("atom"+ " "+i);
			//System.out.println("implicit"+ " "+mol.getAtom(i).getImplicitHydrogenCount());
			//System.out.println("connectedbond"+ " "+ordsum(mol,i));
			//System.out.println("valence"+ " "+valences.get(mol.getAtom(i).getSymbol()));
			return false;
		}else{
			return true;
		}
	}
	
	// Returns how many open positions the atom has. 
	public static int opencounter(IAtomContainer mol, int i)throws CloneNotSupportedException, CDKException, IOException{
		int open = valences.get(mol.getAtom(i).getSymbol()).intValue()- ordsum(mol,i) - mol.getAtom(i).getImplicitHydrogenCount(); 
		return open;
	}

	//Returns the extension bond list of an atom used for class extension.
	//Explanation: Creating bond and checking every step is time consuming and
	//also need to check if there is no bond between pair. Thus, I am just 
	//considering the saturation of the atoms.
	public static  ArrayList<int[]> classext(IAtomContainer mol, int c) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		ArrayList<int[]> bonds = new ArrayList<int[]>();	
		int l=0; //Counting the interactions for the valence check
		if(satcheck(mol,c)){
			for (int i = 0; i < mol.getAtomCount(); i++){
				if(satcheck(mol,i) && c!=i){ //Avoid the self-interactions.
					bonds.add(new int[]{c, i});
					l=l+1;
					if(l==Functions.opencounter(mol, c)){ //Untill the open max are filled.
						break;
					}
				}
			}
		}
		return bonds;
	}
	
	//Different then bond list, find the next vertex.
	public static  int classnext(IAtomContainer mol, int c) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		int target=0;
		for(int i = 0; i < mol.getAtomCount(); i++){
			if(satcheck(mol,c) && satcheck(mol,i) && c!=i){ //Avoid the loops.
				target=target+i;
				break;
			}
		}
		return target;
	}
	
	public static boolean contains(int[] c, int i){
		boolean check= false;
		for(int j=0;j<c.length;j++){
			if(c[j]==i){
				check=true;
			}
		}
		
		return check;
	}
	
	public static boolean containsrec(List<Integer> c, int i){
		boolean check= false;
		for(int j=0;j<c.size();j++){
			if(c.get(j)==i){
				check=true;
			}
		}
		
		return check;
	}
	
	public static int indatom(IAtomContainer mol, IAtom atom){
		int index=0;
		for(int i=0; i<mol.getAtomCount();i++){
			if(mol.getAtom(i)==atom){
				index=index+i;
			}
		}
		return index;
	}
	
	//This function add bonds in the atom saturation based molecular extension.
	public static void conbondadd(IAtomContainer mol, int[] c, int i , int j, Order ord, Map<Integer, Integer> openvls) {
		mol.getAtom(i).setProperty("visit", (int)1);
		mol.addBond(c[j], i, ord);	
		openvls.replace(j, openvls.get(j)-ord.numeric());
		if(contains(c, i)){ 
			openvls.replace(i, openvls.get(i)-ord.numeric());
		}
	}
	
	//This function add bonds in the atom saturation based molecular extension.
	public static void conbondadd1(IAtomContainer mol, List<Integer> c, int i , int j, Order ord, Map<Integer, Integer> openvls) {
		mol.getAtom(i).setProperty("visit", (int)1);
		mol.addBond(c.get(j), i, ord);	
		openvls.replace(j, openvls.get(j)-ord.numeric());
		if(containsrec(c, i)){ 
			System.out.println("contain icinde"+" "+i+" "+openvls.get(i));
			if(openvls.get(i)!=null){
				openvls.replace(i, openvls.get(i)-ord.numeric());
			}
		}
	}
	
	
	/**Similar to previous version; this function add bonds in the atom saturation based molecular extension but
	 * it updates the visit prop in another acontainer. So the atoms can be visited multiple times until the atom saturated. 
	 */
	public static void conbondadd2(IAtomContainer mol, IAtomContainer mol1,int[] c, int i , int j, Order ord, Map<Integer, Integer> openvls) {
		mol1.getAtom(i).setProperty("visit", (int)1);
		mol.addBond(c[j], i, ord);	
		openvls.replace(j, openvls.get(j)-ord.numeric());
		if(contains(c, i)){ 
			openvls.replace(i, openvls.get(i)-ord.numeric());
		}
	}
	
	public static void conbondaddrec(IAtomContainer mol, IAtomContainer mol1,List<Integer> c, int i , int j, Order ord, Map<Integer, Integer> openvls) {
		mol1.getAtom(i).setProperty("visit", (int)1);
		mol.addBond(c.get(j), i, ord);	
		openvls.replace(j, openvls.get(j)-ord.numeric());
		if(containsrec(c, i)){ 
			if(openvls.get(i)!=null){
				openvls.replace(i, openvls.get(i)-ord.numeric());
			}
		}
	}
	//Check the uniqueness of the molecules and store the unique ones in the given hashset
	public static void checkaddupdate(IAtomContainer mol,HashSet<List<Object>> set,List<IAtomContainer> acclass2, Map<Integer, Integer> openvls, Map<Integer, Integer> openvlsfirst) throws CloneNotSupportedException, CDKException, IOException{
		//System.out.println("rep"+ " "+classrep2(mol));
		if(!set.contains(Functions.classrep2(mol))){
			acclass2.add(mol);
			set.add(Functions.classrep2(mol));				
		}
		openvls=openvlsfirst;
	}
	
	//By taking the original mol, the visit property is set to the new atomcontainer used for the recursive structure generation
	public static void setvisitprop(IAtomContainer mol, IAtomContainer acon1){
		for(IAtom atom: mol.atoms()){
			acon1.addAtom(atom);
		}
		
		for(int s=0;s<mol.getAtomCount();s++){
			if(mol.getAtom(s).getProperty("visit")!=null){
				acon1.getAtom(s).setProperty("visit", (int)1);
			}
		}
	}
	
	//Cloning atomcontainers without their bonds but just cloning atoms with their properties. 
	//But it returned some errors like not interacting two atoms with the same atom. Only single interactions.
	public static IAtomContainer containerclone(IAtomContainer mol1) throws CloneNotSupportedException{
		mol1.removeAllBonds();
		IAtomContainer mol2 = mol1.clone();
		//IAtomContainer mol2 = new org.openscience.cdk.AtomContainer();
	return mol2;
	}
	
	/** A new atomcontainer was created but need to have the clone and does not make sense to have a new one and copy the atoms into it.
	 * Rather than that we can simply clone the atomcontainer first but need to create a new atomcontainer not consisting any bonds.
	 */
	public static List<IAtomContainer> acclass= new ArrayList<IAtomContainer>();
	public static HashSet<List<Object>> set=new HashSet<List<Object>>(); 
	public static  List<IAtomContainer> classextcontainer1(IAtomContainer mol, int[] c) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		saturation = new SaturationChecker();
		IAtomContainer molfirst=Functions.containerclone(mol);		
		Map<Integer, Integer> openvlsfirst= Functions.openlist(mol, c);
		Map<Integer, Integer> openvls= Functions.openlist(mol, c);

		if(!saturation.allSaturated(mol)){
			for(int j=0; j<c.length;j++){ //For all the selected atom indices.
				for (int i=0; i<mol.getAtomCount();i++){ // Visit all the atoms in the acontainer
					Order ord= Functions.ordselect(Functions.opendif(openvls.get(j), Functions.opencounter(mol, i)));
					if(ord!=null){
						if(c[j]==i){continue;}
						else{
							if(mol.getAtom(i).getProperty("visit")==null && satcheck(mol,c[j]) && satcheck(mol,i)){ 
								conbondadd2(mol, molfirst,c,i,j,ord,openvls);
								if(visitcheckcount(molfirst)<=molfirst.getAtomCount() && openlistcheck(openvls)){
									//checkaddupdate(mol,set,acclass, openvls, openvlsfirst);
									acclass.add(mol);
									openvls=openvlsfirst;
									classextcontainer1(molfirst,c);					
								}
							}
						}
					}
				}
			}
		}
		return acclass;
	}
	
	public static  List<IAtomContainer> classextrec(IAtomContainer mol, List<Integer> c) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		saturation = new SaturationChecker();
		IAtomContainer molfirst=Functions.containerclone(mol);		
		Map<Integer, Integer> openvlsfirst= Functions.openlistrec(mol, c);
		Map<Integer, Integer> openvls= Functions.openlistrec(mol, c);
		if(!saturation.allSaturated(mol)){
			for(int j=0; j<c.size();j++){ //For all the selected atom indices.
				for (int i=0; i<mol.getAtomCount();i++){ // Visit all the atoms in the acontainer
					Order ord= Functions.ordselect(Functions.opendif(openvls.get(j), Functions.opencounter(mol, i)));
					if(ord!=null){
						if(c.get(j)==i){molfirst.getAtom(i).setProperty("visit", (int)1);}
						else{
							if(mol.getAtom(i).getProperty("visit")==null ){ 
								if(satcheck(mol,c.get(j)) && satcheck(mol,i)){
								conbondaddrec(mol, molfirst,c,i,j,ord,openvls);
								if(visitcheckcount(molfirst)<=molfirst.getAtomCount() && openlistcheck(openvls)){
									//checkaddupdate(mol,set,acclass, openvls, openvlsfirst);
									acclass.add(mol);
									openvls=openvlsfirst;
									classextrec(molfirst,c);					
								}
							}
						}
						}
					}
				}
			}
		}
		return acclass;
	}
	public static HashSet<List<Object>> setall=new HashSet<List<Object>>(); 
	public static  void single(IAtomContainer mol,List<IAtomContainer> acs) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		Functions.propclean(mol);
		System.out.println("represent"+ " "+Functions.cleanatom2impclass2(mol));
		for(String key: Functions.cleanatom2impclass2(mol)){
			//Functions.propclean(mol);
			//acs.addAll(Functions.merged(mol, Functions.atom2impclass(mol).get(key)));
			for(IAtomContainer ac:Functions.merged(mol, Functions.atom2impclass2(mol).get(key))){
				//if(!setall.contains(Functions.classrep2(ac))){
					//setall.add(Functions.classrep2(ac));
					acs.add(ac);
				//}
			}
			System.out.println("singleic"+ " "+acs.size());
		}
	}
	
	public static  List<IAtomContainer>  firstep(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		List<IAtomContainer> acs= new ArrayList<IAtomContainer>();
		Functions.propclean(mol);
		for(int i=0; i<Functions.cleanatom2impclass2(mol).size();i++){
			System.out.println("denetle"+ " "+Functions.atom2impclass2(mol).get(Functions.cleanatom2impclass2(mol).get(i)));
			//Functions.merged(mol,Functions.atom2impclass2(mol).get(Functions.cleanatom2impclass2(mol).get(0)));
			Functions.propclean(mol);
			//System.out.println("sizeabak"+ " "+Functions.merged(mol,Functions.atom2impclass2(mol).get(Functions.cleanatom2impclass2(mol).get(i))).size());
			for(IAtomContainer ac: Functions.merged(mol,Functions.atom2impclass2(mol).get(Functions.cleanatom2impclass2(mol).get(i)))){
				System.out.println("incer");
				acs.add(ac);
				mol.removeAllBonds();
			}
		}
		//System.out.println(Functions.atom2impclass2(mol).get(Functions.cleanatom2impclass2(mol).get(0)));
		
		return acs;
	}
	
	public static List<IAtomContainer> all(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException {
		List<IAtomContainer> acs=new ArrayList<IAtomContainer>();
		List<IAtomContainer> first=firstep(mol);
		acs.addAll(first);
		Functions.propclean(mol);
		//if(visitcheck(mol)){
			for(IAtomContainer ac:first){
				//Functions.propclean(ac);
				single(ac,acs);
				System.out.println("allsizela"+ " "+acs.size());
			}
		//}
		return acs;
	}
	public static List<IAtomContainer> acclass2= new ArrayList<IAtomContainer>();
	public static HashSet<List<Object>> set2=new HashSet<List<Object>>(); 
	//For ignoring the self-interactions.
	public static  List<IAtomContainer> classextcontainer2(IAtomContainer mol, int[] c) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		saturation = new SaturationChecker();
		Map<Integer, Integer> openvlsfirst= Functions.openlist(mol, c);
		Map<Integer, Integer> openvls= Functions.openlist(mol, c);		
		if(!saturation.allSaturated(mol)){
			for(int j=0; j<c.length;j++){ //For all the selected atom indices.
				for(int i=0; i<mol.getAtomCount();i++){ // Visit all the atoms in the acontainer
					Order ord= Functions.ordselect(Functions.opendif(openvls.get(j), Functions.opencounter(mol, i)));
					if(ord!=null){
						if(c[j]==i){ mol.getAtom(i).setProperty("visit", (int)1);}
						else{
							if(mol.getAtom(i).getProperty("visit")==null && satcheck(mol,c[j]) && satcheck(mol,i)){ 
								conbondadd(mol,c,i,j,ord,openvls);
								if(visitcheckcount(mol)<=mol.getAtomCount() && openlistcheck(openvls)){
									IAtomContainer acon1 = new org.openscience.cdk.AtomContainer();
									//checkaddupdate(mol,set2,acclass2,openvls, openvlsfirst);
									acclass2.add(mol);
									openvls=openvlsfirst;
									setvisitprop(mol,acon1);
									classextcontainer2(acon1,c);					
								}
							}
						}
					}
				}
			}
		}	
		return acclass2;
	}
	
	public static  List<IAtomContainer> classextrec2(IAtomContainer mol,List<Integer> c) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		saturation = new SaturationChecker();
		Map<Integer, Integer> openvlsfirst= Functions.openlistrec(mol, c);
		Map<Integer, Integer> openvls= Functions.openlistrec(mol, c);		
		if(!saturation.allSaturated(mol)){
			for(int j=0; j<c.size();j++){ //For all the selected atom indices.
				for(int i=0; i<mol.getAtomCount();i++){ // Visit all the atoms in the acontainer
					Order ord= Functions.ordselect(Functions.opendif(openvls.get(j), Functions.opencounter(mol, i)));
					if(ord!=null){
						if(c.get(j)==i){ mol.getAtom(i).setProperty("visit", (int)1);}
						else{
							if(mol.getAtom(i).getProperty("visit")==null && satcheck(mol,c.get(j)) && satcheck(mol,i)){ 
								conbondadd1(mol,c,i,j,ord,openvls);
								if(visitcheckcount(mol)<=mol.getAtomCount() && openlistcheck(openvls)){
									IAtomContainer acon1 = new org.openscience.cdk.AtomContainer();
									//checkaddupdate(mol,set2,acclass2,openvls, openvlsfirst);
									//if(!set.contains(Functions.classrep2(mol))){
									acclass2.add(mol);				
									//}
									openvls=openvlsfirst;
									setvisitprop(mol,acon1);
									classextrec2(acon1,c);					
								}
							}
						}
					}
				}
			}
		}	
		return acclass2;
	}
	public static List<IAtomContainer> acclass3= new ArrayList<IAtomContainer>();
	public static HashSet<List<Object>> set3=new HashSet<List<Object>>(); 
	//Extending atomcontainer by jumping one index.
	public static  List<IAtomContainer> classextcontainer2jump(IAtomContainer mol, int[] c) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		saturation = new SaturationChecker();
		Map<Integer, Integer> openvlsfirst= Functions.openlist(mol, c);
		Map<Integer, Integer> openvls= Functions.openlist(mol, c);			
		if(!saturation.allSaturated(mol)){
			for(int j=0; j<c.length;j++){ //For all the selected atom indices.
				for (int i=0; i<mol.getAtomCount();i=i+2){ // Visit all the atoms in the acontainer by jumping one index
					Order ord= Functions.ordselect(Functions.opendif(openvls.get(j), Functions.opencounter(mol, i)));
					if(ord!=null){
						if(c[j]==i) {mol.getAtom(i).setProperty("visit", (int)1);}
						else{
							if(mol.getAtom(i).getProperty("visit")==null && satcheck(mol,c[j]) && satcheck(mol,i)){ 
								conbondadd(mol,c,i,j,ord,openvls);
								if(visitcheckcount(mol)<=mol.getAtomCount() && openlistcheck(openvls)){
									IAtomContainer acon1 = new org.openscience.cdk.AtomContainer();
									//checkaddupdate(mol,set3,acclass3,openvls,openvlsfirst);
									acclass3.add(mol);
									openvls=openvlsfirst;
									setvisitprop(mol,acon1);				
									classextcontainer2jump(acon1,c);				
								}
							}
						}
					}
				}
			}
		}	
		return acclass3;
	}
	public static  List<IAtomContainer> classextrec2jump(IAtomContainer mol, List<Integer> c) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		saturation = new SaturationChecker();
		Map<Integer, Integer> openvlsfirst= Functions.openlistrec(mol, c);
		Map<Integer, Integer> openvls= Functions.openlistrec(mol, c);			
		if(!saturation.allSaturated(mol)){
			for(int j=0; j<c.size();j++){ //For all the selected atom indices.
				for (int i=0; i<mol.getAtomCount();i=i+2){ // Visit all the atoms in the acontainer by jumping one index
					Order ord= Functions.ordselect(Functions.opendif(openvls.get(j), Functions.opencounter(mol, i)));
					if(ord!=null){
						if(c.get(j)==i) {mol.getAtom(i).setProperty("visit", (int)1);}
						else{
							if(mol.getAtom(i).getProperty("visit")==null && satcheck(mol,c.get(j)) && satcheck(mol,i)){ 
								conbondadd1(mol,c,i,j,ord,openvls);
								if(visitcheckcount(mol)<=mol.getAtomCount() && openlistcheck(openvls)){
									IAtomContainer acon1 = new org.openscience.cdk.AtomContainer();
									//checkaddupdate(mol,set3,acclass3,openvls,openvlsfirst);
									acclass3.add(mol);
									openvls=openvlsfirst;
									setvisitprop(mol,acon1);				
									classextrec2jump(acon1,c);				
								}
							}
						}
					}
				}
			}
		}	
		return acclass3;
	}
	
	//Make the amount of index jump a parameter of the function
	//Extending atomcontainer by jumping one index.
	public static  List<IAtomContainer> classextjumppar(IAtomContainer mol, int[] c, int next) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		saturation = new SaturationChecker();
		Map<Integer, Integer> openvlsfirst= Functions.openlist(mol, c);
		Map<Integer, Integer> openvls= Functions.openlist(mol, c);			
		if(!saturation.allSaturated(mol)){
			for(int j=0; j<c.length;j++){ //For all the selected atom indices.
				for (int i=0; i<mol.getAtomCount();i=i+next){ // Visit all the atoms in the acontainer by jumping one index
					Order ord= Functions.ordselect(Functions.opendif(openvls.get(j), Functions.opencounter(mol, i)));
					if(ord!=null){
						if(c[j]==i) {mol.getAtom(i).setProperty("visit", (int)1);}
						else{
							if(mol.getAtom(i).getProperty("visit")==null && satcheck(mol,c[j]) && satcheck(mol,i)){ 
								conbondadd(mol,c,i,j,ord,openvls);
								if(visitcheckcount(mol)<=mol.getAtomCount() && openlistcheck(openvls)){
									IAtomContainer acon1 = new org.openscience.cdk.AtomContainer();
									//checkaddupdate(mol,setglo,acclass,openvls,openvlsfirst);
									setvisitprop(mol,acon1);				
									classextjumppar(acon1,c, next);				
								}
							}
						}
					}
				}
			}
		}	
		return acclass;
	}
		
	//Merging two atomcontainer lists
	public static List<IAtomContainer> aclistsum(List<IAtomContainer> acs1, List<IAtomContainer> acs2){
		for(IAtomContainer ac: acs1){
			if(!acs2.contains(ac)){
				acs2.add(ac);
			}
		}
		return acs2;
	}
	
	//Cleaning properties of atoms.
	public static IAtomContainer propclean(IAtomContainer mol){
		for(IAtom atom: mol.atoms()){
			if(atom.getProperty("visit")!=null){
				atom.setProperty("visit", null);
			}
		}
		//mol.removeAllBonds();
		return mol;
	}
	
	//In the merged form of different cases, it cleanes the duplicates.
	public static HashSet<List<Object>> setfinal=new HashSet<List<Object>>();
	public static List<IAtomContainer> cleanduplicates(List<IAtomContainer> ac1) throws CloneNotSupportedException, CDKException, IOException{
		List<IAtomContainer> acclean= new ArrayList<IAtomContainer>();
		for(IAtomContainer ac: ac1){
			if(!setfinal.contains(Functions.classrep(ac))){
				setfinal.add(Functions.classrep(ac));
				acclean.add(ac);
			}
		}
		return acclean;
	}
	
	public static List<IAtomContainer> merged(IAtomContainer mol, List<Integer> c) throws CloneNotSupportedException, CDKException, IOException {
		List<IAtomContainer> total= new ArrayList<IAtomContainer>();
		propclean(mol);
		//mol.removeAllBonds();
		List<IAtomContainer> ac1=Functions.classextrec(mol, c);
		//int olc=0;		
        for(IAtomContainer ac: ac1){
        	//if(Functions.subacsatur(ac)==true){
        		total.add(ac);
            	//Functions.depict(ac, "C:\\Users\\MAY\\Desktop\\aziz"+olc+".png");
            	//olc=olc+1;
        	//}
        } 
		propclean(mol);
		IAtomContainer mol2=mol.clone();
		mol2.removeAllBonds();
		//mol.removeAllBonds();
		List<IAtomContainer> ac2=Functions.classextrec2(mol2, c);
		//int olcu=0;		
        for(IAtomContainer acs: ac2){
        	//if(Functions.subacsatur(ac)==true){
        		total.add(acs);
        		//Functions.depict(acs, "C:\\Users\\MAY\\Desktop\\elis"+olcu+".png");
        		//olcu=olcu+1;
        	//}
        } 
		propclean(mol);
		IAtomContainer mol3=mol.clone();
		mol3.removeAllBonds();
		//mol.removeAllBonds();
		List<IAtomContainer> ac3=Functions.classextrec2jump(mol3, c);
		//int olcum=0;		
        for(IAtomContainer acl: ac3){
        	//if(Functions.subacsatur(ac)==true){
        		total.add(acl);
        		//Functions.depict(acl, "C:\\Users\\MAY\\Desktop\\brosis"+olcum+".png");
        		//olcum=olcum+1;
        	//}
        } 
		//return cleanduplicates(total);
        return total;
	}
	
	//Check whether there is non-visited atom.
	public static boolean visitcheck(IAtomContainer mol) throws CDKException{
		boolean check=false;
		for(int i=0; i<mol.getAtomCount(); i++){
			if(mol.getAtom(i).getProperty("visit")==null){
				check=true;
			}
		}
		return check;
	}
	
	//Set the atomvisit of the given indices.
	public static Map<Integer, Integer> visitatomset(IAtomContainer mol, int[] c) throws CDKException{
		Map<Integer, Integer> map= new HashMap<Integer, Integer>();
		for(int i=0; i< mol.getAtomCount();i++){
			map.put(i, 0);
		}
		return map;
	}
	
	//Check whether there is non-visited atom.
	public static int visitcheckcount(IAtomContainer mol) throws CDKException{
		int count=0;
		for(int i=0; i<mol.getAtomCount(); i++){
			if(mol.getAtom(i).getProperty("visit")!=null){
				count=count+1;
			}
		}
		return count;
	}
		
	//Check whether there is no more open sites.
	public static boolean openlistcheck(Map<Integer, Integer> openvls) throws CDKException{
		boolean check=true;
		for(Integer key: openvls.keySet()){
			if(openvls.get(key)!=0){
				//System.out.println("keys"+ " "+openvls.get(key));
				check=false;
			}
		}
		return check;
	}
	
	//Takes the indices of the atoms to calculate their open values.
	public static Map<Integer, Integer> openlist(IAtomContainer mol, int[] c) throws CloneNotSupportedException, CDKException, IOException{
		Map<Integer, Integer> map= new HashMap<Integer, Integer>();
		for(int i=0; i<c.length;i++){
			map.put(i,Functions.opencounter(mol, c[i]));
		}
		return map;
	}
	
	//Takes the indices of the atoms to calculate their open values.
	public static Map<Integer, Integer> openlistrec(IAtomContainer mol, List<Integer> c) throws CloneNotSupportedException, CDKException, IOException{
		Map<Integer, Integer> map= new HashMap<Integer, Integer>();
		for(int i=0; i<c.size();i++){
			map.put(i,Functions.opencounter(mol, c.get(i)));
		}
		return map;
	}
		
	//If needed, the first values can be also kept in the map for comparing the initial and current values during the extension.
	
	//Takes the indices of the atoms to calculate their open values.
	public static Map<Integer, Integer[]> openlisttotal(IAtomContainer mol, int[] c) throws CloneNotSupportedException, CDKException, IOException{
		Map<Integer, Integer[]> map= new HashMap<Integer, Integer[]>();
		for(int i=0; i<c.length;i++){
			map.put(i,new Integer [] {Functions.opencounter(mol, c[i]),Functions.opencounter(mol, c[i])});
		}
		return map;
	}
		
	
	//To get the difference  between the open sites.
	public static int opendif( int i, int j) throws CDKException {
		if(i>j){
			return j;
		}else{
			return i;
		}
	}
	
	//To choose the order size based on the open differences.
	public static Order ordselect( int i) throws CDKException{
		Order ord=null;
		if(i==1){
			ord= Order.SINGLE;
		}else if(i==2){
			ord= Order.DOUBLE;
		}else if(i>=3){
			ord= Order.TRIPLE;
		}
		return ord;
	}
	
	//To remove the newly added bonds to go back to the original acontainer.
	public static void oriac( IAtomContainer mol) throws CDKException{
		for(IBond bond: mol.bonds()){
			if((int)bond.getProperty("added")==1){
				//System.out.println(bond.getOrder().numeric());
				mol.removeBond(bond);
			}
		}
	}
	
	//Increase bond order for vertex paired bonds int[].
	public static void increaseorder(IAtomContainer mol, int[] bond)throws CloneNotSupportedException, CDKException, IOException {
		IBond add = mol.getBond(mol.getAtom(bond[0]), mol.getAtom(bond[1])); //bond is a 1D array with two entries; edge vertices
		if(add == null){ //Bondmanipulator returns nullpointerexception.					
			mol.addBond(bond[0], bond[1], IBond.Order.SINGLE);
		}
		else{
			BondManipulator.increaseBondOrder(add); 
		}
	}
	
	//Increase bond order of a given pair of edge vertices.
	public static void increaseorder2(IAtomContainer mol, int f, int l)throws CloneNotSupportedException, CDKException, IOException {
		IBond add = mol.getBond(mol.getAtom(f), mol.getAtom(l)); //bond is a 1D array with two entries; edge vertices
		if(add == null){ //Bondmanipulator returns nullpointerexception.					
			mol.addBond(f, l, IBond.Order.SINGLE);
		}
		else{
			BondManipulator.increaseBondOrder(add); 
		}
	}
	
	//Increase bond order of a given pair of edge vertices without using increaseBondOrder.
	public static void increaseorder3(IAtomContainer mol, int f, int l)throws CloneNotSupportedException, CDKException, IOException {
		IBond add = mol.getBond(mol.getAtom(f), mol.getAtom(l)); //bond is a 1D array with two entries; edge vertices
		if(add == null){ //Bondmanipulator returns nullpointerexception.					
			mol.addBond(f, l, IBond.Order.SINGLE);
		}
		else{
		Order order = add.getOrder();	
			if(order == IBond.Order.SINGLE){
				mol.getBond(mol.getAtom(f), mol.getAtom(l)).setOrder(IBond.Order.DOUBLE);
			}
			else if(order == IBond.Order.DOUBLE){
				mol.getBond(mol.getAtom(f), mol.getAtom(l)).setOrder(IBond.Order.TRIPLE);
			}
		}
	}
	
	//Adding bonds to the molecule based on a list of bonds for extension.
	public static IAtomContainer bondadder(IAtomContainer mol, int c)throws CloneNotSupportedException, CDKException, IOException {
		ArrayList<int[]> ext = Functions.classext(mol,c);
		for(int[] bond : ext){
			Functions.increaseorder(mol, bond);
		}
		return mol;
	}	
	
	//Adding bonds to the molecule based on a target atom for extension.
	public static IAtomContainer bondadder2(IAtomContainer mol, int c)throws CloneNotSupportedException, CDKException, IOException {
		int target = Functions.classnext(mol,c);
		if(Functions.satcheck(mol, c)){
			Functions.increaseorder3(mol, c, target);
		}
	return mol;
	}	
	
	//Checks whether the atomcontainer consists a saturated sub atomcontainer.
	public static boolean subacsatur (IAtomContainer acontainer) throws CloneNotSupportedException, CDKException, IOException {
        boolean check = true;
		for(IAtomContainer mol:ConnectivityChecker.partitionIntoMolecules(acontainer).atomContainers()){
        	saturation = new SaturationChecker();
    		if(saturation.isSaturated(mol)==true){
    			check =false;
    		}
        }
		return check;
	}
	
	//This is with the CNI check as explained in the paper.
	public static void gen(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		ArrayList<Long> keys=Functions.cleansymlistkey(mol);
		for(Long key: keys){ //For all the eclasses
			if(intinclass(mol,key) && intwithtarget(mol,key)){ //CNI check. Please look at the paper or the explanation in the test class.
				Functions.gens(mol);
			}
			else{
				Functions.genl(mol);
			}
		}
	}
		
	//Generate graphs by adding new edges to all the elements of the eclasses.
	public static void genl(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		saturation = new SaturationChecker();
		ListMultimap<Long, Long> map= Functions.symlist(mol); //From the atomcontainer to symclass
		for(Long key: Functions.cleansymlistkey(mol)){ //For all the eclasses //cleaning the key reputations in the MultiMap
			for(Long val:map.get(key)){ //For all the values of the classes, which are atom indices.
				Functions.bondadder2(mol, val.intValue()-1); //Atom index out of bounds: 0 <= 8 < 8 error received sym class starts with 1.
				if(!saturation.isSaturated(mol) && subacsatur(mol)){ //Checks whether there is a saturated subgraph or not.
					genl(mol);
				}else{
					//System.out.println("saturated");
				}
			}
		}
	}
	
	//Generate graphs by adding new edges to the target elements of the class representatives.
	public static void gens(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		saturation = new SaturationChecker();
		if(saturation.isSaturated(mol)){
			//System.out.println("saturated");
		}else{
			ListMultimap<String, Integer> map= Functions.atom2symind(mol); //From the atomcontainer to symclass
			for(String key: map.keys()){ //Just extend from the sym representatives represented by their indices. 
				for(Integer val:map.get(key)){
					Functions.bondadder2(mol, val); //if we will use atom2sym, we have long values then; Atom index out of bounds: 0 <= 8 < 8 error received sym class starts with 1.
					if(subacsatur(mol)){ //Checks whether there is a saturated subgraph or not.
						grouprepresentation.Functions.gens(mol);
					}
				}
			}
		}
	}
	
	//With storing the the atomcontainer at each step as described in the paper; generate structures recursively.
	public static List<IAtomContainer> ac= new ArrayList<IAtomContainer>();
	public static void genstore(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		if(Functions.getInstance().aclist.contains(mol)){
			Functions.getInstance().aclist.remove(mol);
		}
		saturation = new SaturationChecker();
		if(saturation.isSaturated(mol)){
			//System.out.println("saturated");
		}else{
			ListMultimap<String, Integer> map= Functions.atom2symind(mol); //From the atomcontainer to symclass
			for(String key: map.keys()){ //Just extend from the sym class representatives 
				for(Integer val:map.get(key)){
					IAtomContainer mol2=Functions.bondadder2(mol, val); //if we will use atom2sym, we have long values then; Atom index out of bounds: 0 <= 8 < 8 error received sym class starts with 1.
					//if(subacsatur(mol2)){ //Checks whether there is a saturated subgraph or not.
						Functions.getInstance().aclist.add(mol2);
					//}
				}
				for(IAtomContainer a: Functions.getInstance().aclist){
					Functions.getInstance().aclist.remove(a);
					genstore(a);
					
				}
			}
		}
	}
	
/**
 * In progress.	
 * 
 */
	
	public static void equivalentclass( IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException {
		//Outofmemory:Java heap space for MolecularFormulaManipulator
		//Equivalent class partitioning
		EquivalentClassPartitioner ec = new EquivalentClassPartitioner(mol);
	    
		//The atoms are stored with their indices in the molecule
		double[] node= new double[mol.getAtomCount()];
	    
	    for(int k=0; k<mol.getAtomCount();k++){
	    	node[k]=k;
	    }
	    
	    //These matrices are computed just for the getEquivalentClass function
	    
	    double[][] bondmat=ec.buildBondMatrix();
	    double[][] nodemat=ec.buildNodeMatrix(node);
	    double[] weight= ec.buildWeightMatrix(nodemat, bondmat);
	    
	    //why we have 8 entries for 7 nodes ? [7, 1, 2, 3, 4, 5, 6, 7] and what these numbers mean ? [6, 1, 2, 3, 4, 5, 6]
	    int[] equivalent= ec.getEquivalentClass(weight);
	    //int[] classes = ec.getTopoEquivClassbyHuXu(mol);
	    //System.out.println("The equivalent classes are:");
	    //System.out.println(Arrays.toString(equivalent));
	}
	

	public static Map<String, Integer> valences; 
	static {

		//The atom valences from CDK.
		valences = new HashMap<String, Integer>();
		
		valences.put("C", 4);
		valences.put("N", 5);
		valences.put("O", 2);
		valences.put("S", 6);
		valences.put("P", 5);
		valences.put("F", 1);
		valences.put("I", 7);
		valences.put("Cl", 5);
		valences.put("Br", 5);
		valences.put("H", 1);
	}
	
}
