package grouprepresentation;

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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.openscience.cdk.interfaces.IBond.Order;
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

public class Functions {
	private static SaturationChecker saturation;

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
		depict.withSize(500, 500).withAtomColors().withAtomNumbers().depict(mol).writeTo(path); 
		
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
				Functions.sintremove(mol, mol.getAtom(i).getSymbol());
			}
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
	
	//Saturation check by comparing the maximum number of atoms can be attached and the currently attached one.
	public static boolean satcheck(IAtomContainer mol, int i)throws CloneNotSupportedException, CDKException, IOException{
		if (mol.getConnectedBondsCount(mol.getAtom(i)) >= (int)valences.get(mol.getAtom(i).getSymbol())){ //I tried with mol.getConnectedBondsCount(mol.getAtom(i))
			return false;
		}else{
			return true;
		}
	}
	
	// Returns how many open positions the atom has. 
	public static int opencounter(IAtomContainer mol, int i)throws CloneNotSupportedException, CDKException, IOException{
		int open = valences.get(mol.getAtom(i).getSymbol()).intValue()- mol.getConnectedBondsCount(mol.getAtom(i)); 
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
					//genl(mol);
				}else{
					System.out.println("done");
				}
			}
		}
	}
	
	//Generate graphs by adding new edges to the target elements of the class representatives.
	public static void gens(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException{
		saturation = new SaturationChecker();
		if(saturation.isSaturated(mol)){
			System.out.println("done");
		}else{
			ListMultimap<String, Integer> map= Functions.atom2symind(mol); //From the atomcontainer to symclass
			for(String key: map.keys()){ //Just extend from the sym representatives of the 
				for(Integer val:map.get(key)){
					Functions.bondadder2(mol, val); //if we will use atom2sym, we have long values then; Atom index out of bounds: 0 <= 8 < 8 error received sym class starts with 1.
					if(subacsatur(mol)){ //Checks whether there is a saturated subgraph or not.
						//grouprepresentation.Functions.gens(mol);
					}
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
	    System.out.println("The equivalent classes are:");
	    System.out.println(Arrays.toString(equivalent));
	}
	

	private static Map<String, Integer> valences; 
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
