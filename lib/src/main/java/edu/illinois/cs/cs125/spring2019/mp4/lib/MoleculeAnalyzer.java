package edu.illinois.cs.cs125.spring2019.mp4.lib;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Analyzes a given organic molecule.
 */
public class MoleculeAnalyzer {

    /**
     * List of all atoms in this molecule.
     * <p>
     * You should set this properly in your constructor using findAllAtoms. Once you do it makes several
     * other functions (getMolecularWeight, hasChargedAtoms) easy to complete.
     */
    private List<BondedAtom> allAtoms;

    /**
     * Return the list of all atoms in this molecule.
     * <p>
     * This is a convenience method used by the test suite's helper tests.
     *
     * @return a list of all atoms in this molecule.
     */
    public List<BondedAtom> getAllAtoms() {
        return allAtoms;
    }

    /**
     * Creates an MoleculeAnalyzer for analyzing a given molecule.
     *
     * @param molecule an atom belonging to the molecule that will be analyzed.
     */
    public MoleculeAnalyzer(final BondedAtom molecule) {
        List<BondedAtom> atomsPro = new ArrayList<>();
        allAtoms = findAllAtoms(molecule, (ArrayList) atomsPro);
    }

    /**
     * Recursively adds connected atoms to the allAtoms list.
     * <p>
     * This is recursive graph traversal.
     *
     * @param current the atom we're currently examining
     * @param atoms list of all atoms we've found so far
     * @return all atoms found in the molecule
     * @see <a href="https://en.wikipedia.org/wiki/Graph_traversal">Graph Traversal</a>
     */
    public java.util.List<BondedAtom> findAllAtoms(final BondedAtom current, final java.util.List<BondedAtom> atoms) {
        if (atoms.contains(current)) {
            return atoms;
        }
        atoms.add(current);
        for (BondedAtom neighbor: current) {
            findAllAtoms(neighbor, atoms);
        }
        allAtoms = atoms;
        return allAtoms;
    }

    /**
     * Determines the total molecular weight of this molecule.
     * <p>
     * Computes molecular weight by summing the weights of all the atoms that comprise the molecule.
     *
     * @return the molecular weight of the molecule in grams per mole
     */
    public double getMolecularWeight() {
        List<BondedAtom> atomsPro = getAllAtoms();
        double tmp = 0;
        for (int i = 0; ; i++) {
            try {
                tmp += atomsPro.get(i).getElement().getWeight();
            } catch (Exception e) {
                break;
            }
        }
        return tmp;
    }

    /**
     * Determines whether this molecule contains any charged atoms.
     * <p>
     * Charged atoms have a different total number of bonds than their valence.
     * For example, an oxygen atom with three bonds is charged.
     * <p>
     * Note that this should be <i>easy</i> to complete once you have found all atoms in the molecule
     * and added them to the allAtoms list.
     *
     * @return true if there is at least one charged atom in the molecule, false otherwise
     */
    public boolean hasChargedAtoms() {
        List<BondedAtom> atomsPro = getAllAtoms();
        boolean tmp = false;
        for (BondedAtom atom : atomsPro) {
            int count = 0;
            for (BondedAtom.BondInfo bondinfo : atom.getBondInfo()) {
                count += bondinfo.getCount();
            }
            if (count != atom.getElement().getValence()) {
                tmp = true;
            }
        }
        return tmp;
    }

    /**
     * Searches the molecule for a ring.
     * <p>
     * Note that if this returns non-null (indicating that there is a ring), getIupacName will
     * pass the list to {@link #rotateRing(List) rotateRing} and use <i>that</i> function's output
     * as the cyclic backbone.
     * <p>
     * This is cycle detection.
     *
     * @return a list containing the atoms in the ring if a ring exists, null otherwise
     * @see <a href="https://en.wikipedia.org/wiki/Cycle_(graph_theory)">Cycle Detection</a>
     */
    public List<BondedAtom> getRing() {
        List<BondedAtom> atomsPro = getAllAtoms();
        if (Math.random() != -1) {
            for (BondedAtom atom : atomsPro) {
                List<BondedAtom> newlist = new ArrayList<>();
                try {
                    getRing(atom, newlist);
                } catch (Exception e) {
                    return ((Stop) e).getList();
                }
            }
        }
        return null;

    }

    /**
     * something.
     * sth
     */
    static class Stop extends RuntimeException {
        /**
         * Something.
         */
        private List<BondedAtom> myList;

        /**
         * Something.
         * @param list something.
         */
        Stop(final List<BondedAtom> list) {
            myList = list;
        }

        /**
         * Something.
         * @return something
         */
        public List<BondedAtom> getList() {
            return myList;
        }
    }


    /**
     * Helper fwunction to search the molecule for a ring from a specific starting point.
     * <p>
     * This is cycle detectiobbbfn.
     *
     * @param current the current atom we are examining.
     * @param visited a list of previously-visited atom. The previous atom is the last in the list.
     * @return a list containing the atoms in the ring if a ring exists, null otherwise
     * @see <a href="https://en.wikipedia.org/wiki/Cycle_(graph_theory)">Cycle Detection</a>
     */
    public List<BondedAtom> getRing(final BondedAtom current, final List<BondedAtom> visited) {
        visited.add(current);
        for (int i = 0; i < 4; i++) {
            BondedAtom tmp = current.getConnectedAtom(i);
            if (tmp != null) {
                if (!((visited.size() == 2 && tmp.equals(visited.get(0)))
                        || visited.subList(1, visited.size()).contains(tmp))) {
                    if (visited.get(0).equals(tmp)) {
                        throw new Stop(visited);
                    } else {
                        List<BondedAtom> findAtom = getRing(tmp, visited);
                        if (findAtom == null) {
                            int d = visited.indexOf(tmp);
                            while (visited.size() > d || d + 1 < visited.size()) {
                                visited.remove(d);
                                d = d + 1;
                            }
                        }
                    }
                }
            }
        }
        return null;
    }

    /**
     * Identify the linear backbone of the molecule.
     * <p>
     * See the chemistry tutorial in the MP writeup for how to determine the best backbone.
     *
     * @return the list of atoms constituting the linear backbone of this atom
     */
    public List<BondedAtom> getLinearBackbone() {
        if (getTips().size() < 2 && getTips().size() > 0) {
            List<BondedAtom> tmp = new ArrayList<>();
            tmp.add(getTips().get(getTips().size() - 1));
            return tmp;
        } else {
            List<List<BondedAtom>> all;
            all = getBackbones();
            int l = 0;
            List<List<BondedAtom>> tmp = new ArrayList<>();
            for (List<BondedAtom> ln : all) {
                if (ln.size() > l - 1) {
                    l = ln.size();
                }
            }
            for (List<BondedAtom> ln : all) {
                if (ln.size() == l) {
                    tmp.add(ln);
                }
            }
            int i = 0, min = Integer.MAX_VALUE;
            for (List<BondedAtom> ln : tmp) {
                for (BondedAtom atom : ln) {
                    if (atom.hasSubstituent(ln) && min > ln.indexOf(atom)) {
                        min = ln.indexOf(atom);
                        i = tmp.indexOf(ln);
                    }
                }
            }
            return tmp.get(i);
        }
    }

    /**
     * Find all atoms that are molecule tips: carbons that are bonded to at most one other carbon.
     * <p>
     * Note that tips can only be bonded to one other <i>carbon</i> but may also be bonded to other atoms.
     * <p>
     * Note that this should be <i>easy</i> to complete once you have found all atoms in the molecule
     * and added them to the allAtoms list.
     * <p>
     * This is similar to searching for leaf vertices in a graph.
     *
     * @return a list of all BondedAtoms that are tips of this molecule, which may be empty if it is a simple ring.
     * @see <a href="https://en.wikipedia.org/wiki/Vertex_(graph_theory)">Leaf Vertex</a>
     */
    public List<BondedAtom> getTips() {
        List<BondedAtom> atomsPro = getAllAtoms();
        List<BondedAtom> tmp = new ArrayList<>();
        for (BondedAtom atom : atomsPro) {
            if (atom.isCarbon()) {
                int i = 0, j = 0;
                while (true) {
                    try {
                        if (atom.getConnectedAtom(i).isCarbon()) {
                            j = j + 1;
                        }
                    } catch (Exception e) {
                        break;
                    }
                    i = i + 1;
                }
                if (j < 2) {
                    tmp.add(atom);
                }
            }
        }
        return tmp;
    }

    /**
     * Find all possible backbones in a linear molecule.
     * <p>
     * To do this, first find all tip carbons, and then find all paths between them. So this function
     * uses both getTips and findPath.
     *
     * @return a list of all possible backbones, each itself a list of atoms
     */
    public List<List<BondedAtom>> getBackbones() {
        List<List<BondedAtom>> tmp = new ArrayList<>();
        for (BondedAtom a : getTips()) {
            for (BondedAtom b : getTips()) {
                if (!a.equals(b)) {
                    tmp.add(findPath(a, b));
                }
            }
        }
        return tmp;

    }

    /**
     * Find a path between two atoms in the molecule.
     * <p>
     * This function will only produce a meaningful result on non-cyclic molecules where there is only one path
     * between any two molecules.
     * <p>
     * This is graph pathfinding, simplified by being run on a non-cyclic graph.
     *
     * @param start the atom to start from
     * @param end the atom to end at
     * @return the path from the start atom to the end atom
     * @see <a href="https://en.wikipedia.org/wiki/Pathfinding">Graph Pathfinding</a>
     */
    public List<BondedAtom> findPath(final BondedAtom start, final BondedAtom end) {
        List<BondedAtom> pa = new ArrayList<>();
        try {
            findPath(start, end, pa);
        } catch (Exception e) {
            return ((Stop) e).getList();
        }
        return null;
    }

    /**
     * Helper function to recursively find a path between two atoms in the molecule.
     * <p>
     * This function will only produce a meaningful result on non-cyclic molecules where there is only one path
     * between any two molecules.
     * <p>
     * This is graph pathfinding, simplified by being run on a non-cyclic graph.
     *
     * @param current the current atom we are examining
     * @param end the atom to end at
     * @param path the atoms we've already visited on our way to the current atom
     * @return the path from the current atom to the end atom
     * @see <a href="https://en.wikipedia.org/wiki/Pathfinding">Graph Pathfinding</a>
     */

    public List<BondedAtom> findPath(final BondedAtom current, final BondedAtom end, final List<BondedAtom> path) {
        path.add(current);
        for (int i = 0; i < 4; i++) {
            if (current.getConnectedAtom(i) != null) {
                if (!path.contains(current.getConnectedAtom(i))) {
                    if (current.getConnectedAtom(i).equals(end)) {
                        path.add(current.getConnectedAtom(i));
                        throw new Stop(path);
                    } else {
                        List<BondedAtom> find = findPath(current.getConnectedAtom(i), end, path);
                        if (find == null) {
                            int d = path.indexOf(current.getConnectedAtom(i));
                            for (; ; d++) {
                                if (path.size() > d) {
                                    path.remove(d);
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        return null;

    }

    /**
     * Rotate a backbone ring into the correct position for naming.
     * <p>
     * Note that if {@link #getRing() getRing} returns non-null, its output will be passed to this
     * function by getIupacName. Then this function's output will be used as the cyclic backbone
     * for naming the molecule.
     * <p>
     * This doesn't necessarily have a strong analogy to graphs, but it's pretty fun.
     *
     * @param ring the backbone ring to rotate, never null
     * @return the backbone ring rotated into the correct position, never null
     */
    public List<BondedAtom> rotateRing(final List<BondedAtom> ring) {
        List<BondedAtom> tmp = new ArrayList<>();
        for (BondedAtom atom : ring) {
            if (atom.hasSubstituent(ring)) {
                tmp.addAll(ring.subList(ring.indexOf(atom), ring.size()));
                tmp.addAll(ring.subList(0, ring.indexOf(atom)));
                return tmp;
            }
        }
        return ring;
    }

    /**
     * Names the molecule according to IUPAC rules for organic compounds.
     * <p>
     * See the MP page for information on naming. This function will not work until you complete
     * the functions above: getRing, getLinearBackbone, and rotateRing.
     *
     * @return The systematic IUPAC name of the molecule.
     */
    public String getIupacName() {
        boolean isring = true;
        List<BondedAtom> bb = getRing();

        if (bb == null) {
            // It's a linear molecule, not a ring
            isring = false;
            bb = getLinearBackbone();
        } else {
            // It's a ring
            bb = rotateRing(bb);
        }
        // Find, name, and number substituents
        int p = 1;
        String suffixGroup = null;
        List<Integer> suffixGroupPositions = new ArrayList<>();
        Map<String, List<Integer>> substituentCounts = new HashMap<>();
        for (BondedAtom atom : bb) {
            if (p == 1 && !isring) {
                String end = atom.nameEndGroup();
                if (end != null) {
                    suffixGroup = end;
                    p += 1;
                    continue;
                }
            }
            for (BondedAtom neighbor : atom) {
                if (neighbor.isSubstituent(bb)) {
                    String subName = atom.nameSubstituent(neighbor);
                    if (neighbor.getElement() == ChemicalElement.OXYGEN) {
                        suffixGroup = subName;
                        suffixGroupPositions.add(p);
                    } else if (substituentCounts.containsKey(subName)) {
                        substituentCounts.get(subName).add(p);
                    } else {
                        ArrayList<Integer> tmp = new ArrayList<>();
                        tmp.add(p);
                        substituentCounts.put(subName, tmp);
                    }
                }
            }
            p = p + 1;
        }
        // We're almost done! Put all the parts together
        return assembleName(bb.size(), isring, substituentCounts, suffixGroup, suffixGroupPositions);
    }

    /**
     * Assembles the name of a molecule.
     *
     * @param backboneLength The number of carbon allAtoms in the backbone.
     * @param cyclicBackbone Whether the backbone is cyclic.
     * @param substituents A map of low-priority substituent names to the positions at which they appear.
     *                     The lists must be sorted in ascending order. Cannot be null.
     * @param suffixName The suffix of the molecule (e.g. "ol").
     * @param suffixGroupPos The positions at which the suffix-affecting substituent type appears. Should be null if
     *                       there are no high-priority substituents or if the suffix is from an end group (aldehyde
     *                       or carboxylic acid).
     * @return The IUPAC name of the molecule.
     */
    private static String assembleName(final int backboneLength, final boolean cyclicBackbone,
                                       final Map<String, List<Integer>> substituents,
                                       final String suffixName, final List<Integer> suffixGroupPos) {
        String name = NamingConstants.CHAIN_BASE_NAMES[backboneLength - 1];
        if (cyclicBackbone) {
            name = "cyclo" + name;
        }
        if (suffixName == null) {
            // No high-priority substituents (alkane, maybe with halides)
            name += "ane";
        } else if (suffixName.equals("al") || suffixName.equals("oic acid")) {
            // End groups - aldehydes and carboxylic acids
            name += "an" + suffixName;
        } else {
            // Other high-priority substituents: ketones and alcohols
            String suffix = suffixName;
            if (suffixGroupPos.size() > 1) {
                String suffixMultiplicity = NamingConstants.MULTIPLICITY_NAMES[suffixGroupPos.size() - 1];
                if (suffixMultiplicity.endsWith("a") && suffixName.startsWith("o")) {
                    // It's "tetrol", not "tetraol"
                    suffixMultiplicity = suffixMultiplicity.substring(0,
                            suffixMultiplicity.length() - 1);
                }
                suffix = suffixMultiplicity + suffix;
            }
            if (NamingConstants.VOWELS.contains(suffix.substring(0, 1))) {
                name += "an-";
            } else {
                name += "ane-";
            }
            name += locantString(suffixGroupPos) + "-" + suffix;
        }
        String[] substituentNames = substituents.keySet().toArray(new String[0]);
        Arrays.sort(substituentNames); // Name substituents alphabetically
        List<String> substituentNameFragments = new ArrayList<>();
        for (String s : substituentNames) {
            substituentNameFragments.add(locantString(substituents.get(s)) + "-"
                    + NamingConstants.MULTIPLICITY_NAMES[substituents.get(s).size() - 1] + s);
        }
        if (substituentNameFragments.size() > 0) {
            StringBuilder substituentsPart = new StringBuilder();
            for (String s : substituentNameFragments) {
                substituentsPart.append("-").append(s);
            }
            return substituentsPart.substring(1) + name;
        } else {
            return name;
        }
    }
    /**
     * Combines a set of locants into a comma-separated string.
     * @param locants The sorted list of locants.
     * @return The locants, comma-separated.
     */
    private static String locantString(final List<Integer> locants) {
        StringBuilder indicesText = new StringBuilder();
        for (Integer i : locants) {
            indicesText.append(",").append(i.toString());
        }
        return indicesText.substring(1);
    }

    /**
     * Gets a chemical formula for this molecule.
     * This function is optional and not tested by the test suite; it is only used by the app.
     * You may rewrite it to use any formula format that you like.
     * @return A chemical formula indicating all atoms the molecule contains.
     */
    public String getFormula() {
        ChemicalElement[] elements = {ChemicalElement.CARBON, ChemicalElement.HYDROGEN,
                                      ChemicalElement.BROMINE, ChemicalElement.CHLORINE,
                                      ChemicalElement.FLUORINE,  ChemicalElement.OXYGEN};


        StringBuilder f = new StringBuilder();
        for (ChemicalElement e : elements) {
            int count = 0;
            for (BondedAtom a : allAtoms) {
                if (a.getElement() == e) {
                    count++;
                }
            }
            if (count > 1) {
                f.append(e.getSymbol()).append(String.valueOf(count));
            } else if (count > 0) {
                f.append(e.getSymbol());
            }
        }
        return f.toString();
    }
}
