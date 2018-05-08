package grouprepresentation;

/**
 * This is the test implementation of the functions developed for re-implementation of Faulon's paper.
 * 
 * 
 * Reference (R): Faulon, Jean Loup. "On using graph-equivalent classes for the structure
 * elucidation of large molecules." Journal of chemical information and computer sciences 
 * 32.4 (1992): 338-348.
 */

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import com.google.common.collect.ListMultimap;

public class Implement {
	// It was assumed that the carbon atoms have such H atoms based on NMR( like in Faulon)
	public static void main(String[] args) throws IOException, CDKException, CloneNotSupportedException, NullPointerException{
		//Adding the hydrogen atoms
		IAtomContainer acontainer = new org.openscience.cdk.AtomContainer();

        acontainer.addAtom(new Atom("C")); //1
        acontainer.addAtom(new Atom("C")); //2
        acontainer.addAtom(new Atom("C")); //3
        acontainer.addAtom(new Atom("H")); //4
        acontainer.addAtom(new Atom("H")); //5
        acontainer.addAtom(new Atom("H")); //6
        acontainer.addAtom(new Atom("H")); //7
        acontainer.addAtom(new Atom("H")); //8
        acontainer.addAtom(new Atom("H")); //9
              
        //Adding the hydrogen atom bonded bonds.
        //acontainer.addBond(0, 3, IBond.Order.SINGLE);
        //acontainer.addBond(0, 4, IBond.Order.SINGLE);
        acontainer.addBond(1, 5, IBond.Order.SINGLE);
        acontainer.addBond(1, 6, IBond.Order.SINGLE);
        acontainer.addBond(2, 7, IBond.Order.SINGLE);
        acontainer.addBond(2, 8, IBond.Order.SINGLE);
        acontainer.addBond(0, 1, IBond.Order.SINGLE);
        acontainer.addBond(0, 2, IBond.Order.SINGLE);
        acontainer.addBond(1, 2, IBond.Order.SINGLE);
        
        
        /**
         * The basic functions are explained below
         */
        
        
        /**
         * First, check whether the molecule is ideal or not.
         * If not, make it ideal.
         */
        boolean check = Functions.idealcheck(acontainer);
        
        if(check==false){
        	Functions.makeIdeal(acontainer);
        }
        
        
        /**
         * The symmetry values of the atoms can be created by calling the canon.sym function, if needed.
         * Symmetry values are computed by using Canon class of CDK. "canon.symmetry"
         */
        
        long[] sym= Functions.canonsym(acontainer);
        System.out.println(Arrays.toString(sym));
        
        /**
         * For makeIdeal, the nonsymdetect function is used. 
         * The non-symmetric ones are detected. 
         * Then the interactions of these atoms are removed.
         */
        
		List<Integer> indk= Functions.nonsymdetect(sym);
		System.out.println(indk);
		
		
        /**
         * The symmetry lists are created with symlist function. For all the symmetry values,
         * their equivalent atoms' indices are listed.
         */
        
        ListMultimap<Long, Long> symlist= Functions.symlist(acontainer);
        System.out.println(symlist);
        
        /**
         * atom2symc functions returns the symmetric class values of elements
         */
        
        ListMultimap<String, Long> atom2sym= Functions.atom2symc(acontainer);
        System.out.println(atom2sym);
        
        /**
         * atom2symind functions returns the symmetric class of elements with their extension indices.
         */
        
        ListMultimap<String, Integer> atom2symind= Functions.atom2symind(acontainer);
        System.out.println(atom2symind);
        
        /**
         * For atoms, their possible interactions are detected by classext.
         * For instance, for atom with index 0, its extendsion list is created.
         */
        
        //ArrayList<int[]> classext = Functions.classext(acontainer, 0);

        
        /**
         * For detecting the potential vertex for interaction, classnext can be used.
         */
       
        int classnext = Functions.classnext(acontainer, 0);
        System.out.println(classnext);
        
        /**
         * To check whether the atomcontainer includes a saturated sub atomcontainer,
         * subacsatur function is coded. It returns boolean value.
         */
        
        boolean subacsatur =Functions.subacsatur(acontainer);
        System.out.println(subacsatur);
        
        /**
         * For generation,different options are available:
         * 
         * 1. genl: Extending the molecule by adding bonds to all symmetric vertices, list by list. List is the list of symmetric vertices.
         * 2. gens: Extending the molecule by adding bonds only to the symmetry class representatives.
         * If vertices are symmetric, one index is taken as the representative.
         * 3. gen : First checks CNI; then decide for generation method.
         * 
         * CNI: Checks isomorphism. If the sym class elements are not interacting with each other and with extension vertex; 
         * then only the sym class repsentatives can be used for extension (gens); if not, then all the symmetric vertices should be 
         * extended (genl).
         */
        
        Functions.gen(acontainer);
	}
}
