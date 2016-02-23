///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines, Joseph   //
// Ramsey, and Clark Glymour.                                                //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA //
///////////////////////////////////////////////////////////////////////////////

package edu.pitt.csb.stability;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import edu.cmu.tetrad.data.ColtDataSet;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.RandomSampler;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.GraphSearch;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.CombinatoricsUtils;

import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.RecursiveTask;

/**
 * Runs a search algorithm over a N subsamples of size b to asses stability
 * as in "Stability Selection" and "Stability Approach to Regularization Selection"
 *
 * This is under construction...likely to be buggy
 *
 * Created by ajsedgewick on 9/4/15.
 */
public class StabilitySearch {

    final private DataSet data;
    final private List<Set<Integer>> samps;
    final private DoubleMatrix2D thetaMat;
    final private HashMap<Edge, Double> edgeStab = new HashMap<>();
    final private DataGraphSearch gs;
    private int numVars;
    private int numSamps;
    private int numSubs;
    private int subSize;
    private boolean cpss;
    private String runDir = null;
    private PrintWriter pwLog = null;
    private PrintWriter pwTime = null;

    final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

    /**
     * Constructor using default subsize and numsubs from huge package in R
     *
     * @param data
     * @param gs
     */
    public StabilitySearch(DataSet data, DataGraphSearch gs){
        this.data = data;
        this.gs = gs;
        this.numSamps = data.getNumRows();
        this.numSubs = 20;
        this.numVars = data.getNumColumns();
        this.subSize = (int) Math.floor(.8*numSamps);
        if(numSamps > 144) subSize = (int) Math.floor(10.0*Math.sqrt(numSamps));
        this.cpss = false;

        this.samps = subSampleNoReplacement(false);//numSamps, numSubs, subSize);
        this.thetaMat = DoubleFactory2D.dense.make(numVars, numVars, 0.0);
        //this(data, gs, 20, subSize);
    }

    /**
     *
     * @param data
     * @param gs
     * @param N Number of subsamples
     * @param b Size of Subsamples
     */
    public StabilitySearch(DataSet data, DataGraphSearch gs, int N, int b){
        this.data = data;
        this.gs = gs;
        this.numSamps = data.getNumRows();
        this.numSubs = N;
        this.numVars = data.getNumColumns();
        this.subSize = b;
        this.cpss = false;

        if(b > data.getNumRows()) {throw new IllegalArgumentException("Subsample size must be less than number of samples");}

        this.samps = subSampleNoReplacement(false);//numSamps, N, b);
        this.thetaMat = DoubleFactory2D.dense.make(numVars, numVars, 0.0);
    }

    /**
     * When the size of subsamples is omitted, we use the CPSS approach of Shah and Samworth (2013)
     *
     * @param data
     * @param gs
     * @param N Number of subsamples
     */
    public StabilitySearch(DataSet data, DataGraphSearch gs, int N, boolean cpss){
        this.data = data;
        this.gs = gs;
        this.numSamps = data.getNumRows();
        this.numVars = data.getNumColumns();
        this.cpss = cpss;

        if(cpss) {
            this.subSize = (int) Math.floor(numSamps / 2.0);
            this.numSubs = N*2;
        } else {
            this.subSize = (int) Math.floor(.8*numSamps);
            if(numSamps > 144) subSize = (int) Math.floor(10.0*Math.sqrt(numSamps));
            this.numSubs = N;
        }

        this.samps = subSampleNoReplacement(cpss);//numSamps, N/2, subSize);
        this.thetaMat = DoubleFactory2D.dense.make(numVars, numVars, 0.0);
    }

    public StabilitySearch(DataSet data, DataGraphSearch gs, List<Set<Integer>> samps){
        this.data = data;
        this.gs = gs;
        this.numSamps = data.getNumRows();
        this.numSubs = samps.size();
        this.numVars = data.getNumColumns();
        this.subSize = samps.get(0).size();
        this.cpss = false;

        this.samps = samps;
        this.thetaMat = DoubleFactory2D.dense.make(numVars, numVars, 0.0);
    }

    /**
     * sets run dir where logs of run restuls, details and parameters are saved
     * @param dir
     */
    public void setRunDir(String dir){
        runDir = dir;
        File rd = new File(dir);
        if(rd.exists())
            throw new IllegalArgumentException("Directory already exists. Please use different name to avoid overwriting");
        else if(!rd.mkdir())
            throw new IllegalStateException("Directory not created for some reason...");
        else {
            (new File(rd, "subnets")).mkdir();
            try{pwLog = new PrintWriter(new File(runDir, "runlog.txt"));
                pwTime = new PrintWriter(new File(runDir, "subnetTimes.txt"));
            }
            catch (Throwable t) { t.printStackTrace();}
            writeLog();
            writeSubSamps(rd);
        }
    }

    /**
     * save log file to run directory
     */
    public void writeLog(){
        try{
            pwLog.println(toString());
            pwLog.flush();
        } catch (Throwable T) {
            T.printStackTrace();
        }

    }

    public String toString(){
        StringWriter sw = new StringWriter();
        sw.write("Num Observations: " + numSamps + "\n");
        sw.write("Num Variables: " + numVars + "\n");
        sw.write("Num subsamples: " + numSubs + "\n");
        sw.write("Subsample size: " + subSize + "\n");
        sw.write("Running cpss?: " + cpss + "\n");
        sw.write("Class of graph search: " + gs.getClass().getName() + "\n");
        sw.write("Graph search parameters: " + Arrays.toString(gs.searchParams) + "\n");
        return sw.toString();
    }

    /**
     * save a file list of subsample indexes to directory f
     */
    public void writeSubSamps(File f){
        File ss = new File(f, "subsamp_inds.txt");
        try{
            PrintWriter pw = new PrintWriter(ss);
            for(Set s : samps){
                pw.println(s);
            }
            pw.flush();
        } catch (Throwable T) {
            T.printStackTrace();
        }

    }

    public void saveStability(){
        if(runDir == null)
            throw new IllegalArgumentException("runDir not found!");

        //save results
        try {
            printEdgeStability(true, new PrintWriter(new File(runDir, "directed_stab.txt")));
            printEdgeStability(false, new PrintWriter(new File(runDir, "undirected_stab.txt")));

            double meanAdj = 0.0;
            double meanDir = 0.0;
            for(Edge e: edgeStab.keySet()){
                double d = edgeStab.get(e);
                meanAdj += d;
                if(e.getEndpoint1()==Endpoint.TAIL && e.getEndpoint2()==Endpoint.TAIL)
                    meanDir += 2*d;
                else if (e.getEndpoint1()==Endpoint.ARROW && e.getEndpoint2()==Endpoint.ARROW)
                    meanDir += 2*d;
                else
                    meanDir += d;

            }
            meanDir/=2.0;
            pwLog.println("Mean number of adjacences: " + meanAdj);
            pwLog.println("Mean number of directions: " + meanDir);
            pwLog.flush();
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

    //returns an adjacency matrix containing the edgewise instability as defined in Liu et al
    public void searchSerial() {

        long startMain = System.currentTimeMillis();
        for (int s = 0; s < numSubs; s++) {
            DataSet dataSubSamp = data.subsetRows(ArrayUtils.toPrimitive(samps.get(s).toArray(new Integer[0])));
            long start = System.currentTimeMillis();
            Graph g = gs.search(dataSubSamp);
            long end = System.currentTimeMillis();

            addEdges(g);
            if (runDir != null) {
                saveNet(s, g);
                pwTime.write("sn" + s + "\t" + (end-start)/1000.0 + "\n");
            }

            //DoubleMatrix2D curAdj = MixedUtils.skeletonToMatrix(g);
            //thetaMat.assign(curAdj, Functions.plus);
        }
        //thetaMat.assign(Functions.mult(1.0 / numSubs));
        //thetaMat.assign(thetaMat.copy().assign(Functions.minus(1.0)), Functions.mult).assign(Functions.mult(-2.0));
        //return thetaMat;
        updateThetaMat();
        if (runDir != null) {
            pwTime.write("Total\t" + (System.currentTimeMillis()-startMain)/1000.0 + "\n");
            pwTime.flush();
        }

    }

    //returns an adjacency matrix containing the edgewise instability as defined in Liu et al
    public  void searchParallel(){

        //final int numVars = data.getNumColumns();
        //final DoubleMatrix2D thetaMat = DoubleFactory2D.dense.make(numVars, numVars, 0.0);
        //final HashMap<Edge, Double> thetaList = new HashMap<>();

        //final int[][] samps = subSampleNoReplacement(data.getNumRows(), b, N);

        class StabilityAction extends RecursiveAction{
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to){
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }

            //could avoid using syncronized if we keep track of array of mats and add at end, but that needs lots of
            //memory
            /*private void addToMat(DoubleMatrix2D matSum, DoubleMatrix2D curMat){
                matSum.assign(curMat, Functions.plus);
                if(!matSum.equals(Algebra.DEFAULT.transpose(matSum))) {
                    System.out.println("NOT SYMMETRIC!!!:\n" + matSum + "\ncurmat:\n" + curMat);
                    //throw new IllegalStateException("not symmetric");
                }
            }*/

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        DataGraphSearch curGs;
                        DataSet dataSubSamp;
                        synchronized (gs) {
                            int[] inds = ArrayUtils.toPrimitive(samps.get(s).toArray(new Integer[0]));
                            dataSubSamp = MixedUtils.deepCopy((ColtDataSet)data.subsetRows(inds));
                            curGs = gs.copy();
                        }

                        long start = System.currentTimeMillis();
                        Graph g = curGs.search(dataSubSamp);
                        long end = System.currentTimeMillis();

                        synchronized (edgeStab){
                            addEdges(g);
                            if(runDir != null) {
                                saveNet(s, g);
                                pwTime.write("sn" + s + "\t" + (end-start)/1000.0 + "\n");
                                pwTime.flush();
                            }
                        }

                        //DoubleMatrix2D curAdj = MixedUtils.skeletonToMatrix(g); //set weights so that undirected stability works
                        //synchronized (thetaMat) {
                        //    addToMat(thetaMat, curAdj);
                        //}
                    }

                    return;
                } else {
                    List<StabilityAction> tasks = new ArrayList<StabilityAction>();

                    final int mid = (to - from) / 2;

                    tasks.add(new StabilityAction(chunk, from, from + mid));
                    tasks.add(new StabilityAction(chunk, from + mid, to));

                    invokeAll(tasks);

                    return;
                }
            }

        }

        long startMain = System.currentTimeMillis();
        final int chunk = 1;

        pool.invoke(new StabilityAction(chunk, 0, numSubs));

        updateThetaMat();
        if (runDir != null) {
            pwTime.write("Total\t" + (System.currentTimeMillis()-startMain)/1000.0 + "\n");
            pwTime.flush();
        }

        //thetaMat.assign(Functions.mult(1.0 / numSubs));

        //do this elsewhere
        //thetaMat.assign(thetaMat.copy().assign(Functions.minus(1.0)), Functions.mult).assign(Functions.mult(-2.0));
        //return thetaMat;
    }

    /**
     * Add edges to edge stability counts, Converts "<->" to "---"
     * @param g
     */
    private void addEdges(Graph g){
        for(Edge e : g.getEdges()){
            //converts bidirected to undirected
            if(e.getEndpoint1()==Endpoint.ARROW && e.getEndpoint2()==Endpoint.ARROW){
                e.setEndpoint1(Endpoint.TAIL);
                e.setEndpoint2(Endpoint.TAIL);
            }

            if(edgeStab.containsKey(e)){
                edgeStab.put(e, edgeStab.get(e) + 1.0/numSubs);
            } else {
                edgeStab.put(e, 1.0/numSubs);
            }
        }
    }

    private void saveNet(int ind, Graph g){
        GraphUtils.saveGraph(g, new File(runDir + "/subnets/sn" + ind + ".txt"), false);
    }

    private void updateThetaMat(){
        String dlm = "\t";
        DecimalFormat df = new DecimalFormat("#.######");
        HashMap<Edge, Double> hm;
        hm = new HashMap<>();
        for(Edge e : edgeStab.keySet()){
            double d = edgeStab.get(e);
            e = new Edge(e);
            e.setEndpoint1(Endpoint.TAIL);
            e.setEndpoint2(Endpoint.TAIL);
            if(hm.containsKey(e)){
                hm.put(e, hm.get(e) + d);
            } else{
                hm.put(e,d);
            }
        }

        for (Edge e : hm.keySet()) {
            int x = data.getColumn(e.getNode1());
            int y = data.getColumn(e.getNode2());
            double d = hm.get(e);
            thetaMat.set(x,y,d);
            thetaMat.set(y,x,d);
        }
    }

    //needs a symmetric matrix
    //array of averages of instability matrix over [all, cc, cd, dd] edges
    public static double[] totalInstabilityUndir(DoubleMatrix2D xi, List<Node> vars){
        if (vars.size()!= xi.columns() || vars.size()!= xi.rows()) {
            throw new IllegalArgumentException("stability mat must have same number of rows and columns as there are vars");
        }

        Algebra al = new Algebra();
        //DoubleMatrix2D xiu = MGM.upperTri(xi.copy().assign(al.transpose(xi)),1);

        DoubleMatrix2D xiu = xi.copy().assign(xi.copy().assign(Functions.minus(1.0)), Functions.mult).assign(Functions.mult(-2.0));

        double[] D = new double[4];
        int[] discInds = MixedUtils.getDiscreteInds(vars);
        int[] contInds = MixedUtils.getContinuousInds(vars);
        int p = contInds.length;
        int q = discInds.length;
        double temp = MGM.upperTri(xiu.copy(),1).zSum();
        D[0] = temp/ ((p+q-1.0)*(p+q)/2.0);
        temp = MGM.upperTri(xiu.viewSelection(contInds, contInds).copy(), 1).zSum();
        D[1] = temp/(p*(p-1.0)/2.0);
        temp = xiu.viewSelection(contInds, discInds).zSum();
        D[2] = temp/(p*q);
        temp = MGM.upperTri(xiu.viewSelection(discInds, discInds).copy(), 1).zSum();
        D[3] = temp/(q*(q-1.0)/2.0);

        return D;
    }

    //array of averages of instability matrix over [all, cc, cd, dd] edges
    public static double[] totalInstabilityDir(DoubleMatrix2D xi, List<Node> vars){
        if (vars.size()!= xi.columns() || vars.size()!= xi.rows()) {
            throw new IllegalArgumentException("stability mat must have same number of rows and columns as there are vars");
        }

        double[] D = new double[4];
        int[] discInds = MixedUtils.getDiscreteInds(vars);
        int[] contInds = MixedUtils.getContinuousInds(vars);
        int p = contInds.length;
        int q = discInds.length;
        D[0] = xi.zSum()/ ((p+q-1)*(p+q)/2);

        D[1] = xi.viewSelection(contInds, contInds).zSum()/(p*(p-1));
        D[2] = xi.viewSelection(contInds, discInds).zSum()/(p*q);
        D[3] = xi.viewSelection(discInds, discInds).zSum()/(q*(q-1));

        return D;
    }

    public HashMap<Edge,Double> edgeStability(boolean directed){
        HashMap<Edge,Double> hm = new HashMap<>();
        for(int i = 0; i < numVars; i++){
            for(int j = 0; j < numVars; j++){
                double curVal = thetaMat.get(i,j);
                if(curVal != 0) {
                    Edge e = new Edge(data.getVariable(i), data.getVariable(j), Endpoint.TAIL, Endpoint.TAIL);
                    hm.put(e, curVal);
                }
            }
        }
        return hm;
    }

    /**
     * print list of edge stabilities to tab-delimited format
     * @param directed
     * @param ps
     */
    public void printEdgeStabilityMat(boolean directed, PrintStream ps){
        String dlm = "\t";
        HashMap<Edge,Double> hm = edgeStability(directed);
        for(Edge e : hm.keySet()){
            String eType = "";
            if(e.getEndpoint1()==Endpoint.TAIL) {
                if (e.getEndpoint2() == Endpoint.TAIL)
                    eType = "---";
            }
            ps.println(e.getNode1() + dlm + e.getNode2() + dlm + eType + dlm + hm.get(e));
        }
    }

    public DoubleMatrix2D getThetaMat(){
        return thetaMat;
    }

    /**
     * print list of edge stabilities to tab-delimited format, round to 5 decimal places. currently does not handle
     * circle endpoints...
     * @param directed
     * @param pw
     */
    public void printEdgeStability(boolean directed, PrintWriter pw){
        String dlm = "\t";
        DecimalFormat df = new DecimalFormat("#.######");
        HashMap<Edge, Double> hm;
        if(directed) {
            hm = edgeStab;
        } else {
            hm = new HashMap<>();
            for(Edge e : edgeStab.keySet()){
                double d = edgeStab.get(e);
                e = new Edge(e);
                e.setEndpoint1(Endpoint.TAIL);
                e.setEndpoint2(Endpoint.TAIL);
                if(hm.containsKey(e)){
                    hm.put(e, hm.get(e) + d);
                } else{
                    hm.put(e,d);
                }
            }
        }

        for (Edge e : hm.keySet()) {
            String eType = edgeType(e);
            double d = hm.get(e);

            pw.println(e.getNode1() + dlm + e.getNode2() + dlm + eType + dlm + df.format(d));
        }
        pw.flush();
    }

    public static String edgeType(Edge e){
        StringBuilder buf = new StringBuilder();

        Endpoint endptTypeA = e.getEndpoint1();
        Endpoint endptTypeB = e.getEndpoint2();

        if (endptTypeA == Endpoint.TAIL) {
            buf.append("-");
        } else if (endptTypeA == Endpoint.ARROW) {
            buf.append("<");
        } else if (endptTypeA == Endpoint.CIRCLE) {
            buf.append("o");
        }

        buf.append("-");

        if (endptTypeB == Endpoint.TAIL) {
            buf.append("-");
        } else if (endptTypeB == Endpoint.ARROW) {
            buf.append(">");
        } else if (endptTypeB == Endpoint.CIRCLE) {
            buf.append("o");
        }

        return buf.toString();
    }

    //returns an numSub by subSize matrix of subsamples of the sequence 1:sampSize
    public List<Set<Integer>> subSampleNoReplacement(boolean cpss){//(int sampSize, int numSub, int subSize){

        if (subSize < 1) {
            throw new IllegalArgumentException("Subsample size must be > 0.");
        }

        List<Set<Integer>> subSamps = new ArrayList<Set<Integer>>(numSubs);
        double numSubSamples = numSubs;
        if(cpss) {numSubSamples = numSubs/2.0;}
        for(int i = 0; i < numSubSamples; i++) {
            Set<Integer> curSamp;// = new HashSet<Integer>(subSize);

            // In cpss,  complementary pairs of subsets are drawn together
            // they aren't explicit about non-replacement, but its implied...
            if(cpss){
                Set<Integer> curSamp2;
                do {
                    curSamp = subSampleIndices(numSamps, subSize);
                    curSamp2 = getComplement(curSamp, numSamps);
                } while (subSamps.contains(curSamp) || subSamps.contains(curSamp2));
                subSamps.add((i*2), curSamp);
                subSamps.add((i*2)+1, curSamp2);

            } else {
                do {
                    curSamp = subSampleIndices(numSamps, subSize);
                } while (subSamps.contains(curSamp));

                subSamps.add(i, curSamp);
            }
        }

        if(cpss && subSamps.size() != numSubSamples*2)
            throw new IllegalStateException("CPSS subsamples not equal to 2*num of requested subsamps");
        else if(!cpss && subSamps.size() != numSubSamples)
            throw new IllegalStateException("subsamples not equal to num of requested subsamps. requested: " + numSubSamples + "found: " + subSamps.size());

        return subSamps;
    }

    public static Set<Integer> subSampleIndices(int N, int subSize){
        List<Integer> indices = new ArrayList<Integer>(N);
        for (int i = 0; i < N; i++) {
            indices.add(i);
        }

        Collections.shuffle(indices);
        HashSet<Integer> samp = new HashSet<Integer>(subSize);
        for(int i = 0; i < subSize; i++){
            samp.add(indices.get(i));
        }
        return samp;
    }

    /**
     * returns a complementary set of integers ranging from 0 to n-1 to s both of size floor n/2
     * if n is odd, one int is not in either set at random
     * @param s
     * @param n
     * @return
     */
    public static Set<Integer> getComplement(Set<Integer> s, int n){
        if(s.size()!=(int) Math.floor(n/2.0)) {throw new IllegalArgumentException("Set s must be of size floor(n/2) found s: " + s.size() + " n: " + n);}
        Set<Integer> c = new HashSet<Integer>(s.size());

        boolean dropOne = n % 2 > 0;
        int skipI = -1;
        if(dropOne)
            skipI = (new Random()).nextInt((int) Math.ceil(n/2.0));

        for(int i = 0; i < n; i++){
            if(!s.contains(i)) {
                if (dropOne && c.size() == skipI) {
                    dropOne = false;
                    continue;
                } else
                    c.add(i);
            }
        }

        if(s.size() != c.size())
            throw new IllegalStateException("Complement is not same size as original!");

        return c;
    }

    private static void tests1(){
        String fn = "/Users/ajsedgewick/tetrad_mgm_runs/run2/networks/DAG_0_graph.txt";
        Graph trueGraph = GraphUtils.loadGraphTxt(new File(fn));
        DataSet ds = null;
        try {
            ds = MixedUtils.loadData("/Users/ajsedgewick/tetrad_mgm_runs/run2/data/", "DAG_0_data.txt");
        } catch (Throwable t) {t.printStackTrace();}

        double lambda = .1;
        //SearchWrappers.MGMWrapper mgm = new SearchWrappers.MGMWrapper(new double[]{lambda, lambda, lambda});
        SearchWrappers.PcStableWrapper mgm = new SearchWrappers.PcStableWrapper(new double[]{.05});
        long start = System.currentTimeMillis();
        StabilitySearch ss = new StabilitySearch(ds, mgm, 20, 200);
        ss.searchSerial();
        DoubleMatrix2D xi = ss.getThetaMat();
        long end = System.currentTimeMillis();
        System.out.println("Not parallel: " + ((end-start)/1000.0));

        start = System.currentTimeMillis();
        ss = new StabilitySearch(ds, mgm, 20, 200);
        ss.searchParallel();
        DoubleMatrix2D xi2 = ss.getThetaMat();
        end = System.currentTimeMillis();
        System.out.println("Parallel: " + ((end-start)/1000.0));

        System.out.println(xi);
        System.out.println(xi2);
    }

    private static DataSet genTestData(){
        Graph g = GraphConverter.convert("X1-->X2,X3-->X2,X4-->X5");
        //Graph g = GraphConverter.convert("X1-->X2");
        //simple graph pm im gen example

        HashMap<String, Integer> nd = new HashMap<>();
        nd.put("X1", 0);
        nd.put("X2", 5);
        nd.put("X3", 5);
        nd.put("X4", 5);
        nd.put("X5", 5);

        g = MixedUtils.makeMixedGraph(g, nd);

        //GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(g, "Split(-1.5,-1,1,1.5)");
        GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(g, "2");
        try{pm.setParameterExpression("s1", "1");
            pm.setParameterExpression("s5", "1");}
        catch(Throwable t) {t.printStackTrace();}

        //GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(g, "1");
        System.out.println(pm);

        GeneralizedSemIm im = MixedUtils.GaussianCategoricalIm(pm, false);
        System.out.println(im);

        int samps = 1000;
        DataSet ds = im.simulateDataAvoidInfinity(samps, true);
        ds = MixedUtils.makeMixedData(ds, nd);
        return ds;
    }

    private static void tests2(){
       DataSet ds = genTestData();

        double lambda = .05;
        //SearchWrappers.MGMWrapper mgm = new SearchWrappers.MGMWrapper(new double[]{lambda, lambda, lambda});
        SearchWrappers.PcStableWrapper mgm = new SearchWrappers.PcStableWrapper(new double[]{.05});
        long start = System.currentTimeMillis();
        int N = 100;
        //StabilitySearch ss = new StabilitySearch(ds, mgm, N, 500);
        StabilitySearch ss = new StabilitySearch(ds, mgm, N, true);
        //ss.samps = new ArrayList<Set<Integer>>();
        ss.setRunDir("/Users/ajsedgewick/stab_test");
        ss.searchParallel();
        ss.saveStability();
        DoubleMatrix2D xi = ss.getThetaMat();
        long end = System.currentTimeMillis();
        double time1 = ((end-start)/1000.0);

        start = System.currentTimeMillis();
        //ss = new StabilitySearch(ds, mgm, N, 500);
        ss = new StabilitySearch(ds, mgm, N, true);
        ss.searchSerial();
        DoubleMatrix2D xi2 =  ss.getThetaMat();
        end = System.currentTimeMillis();

        System.out.println("Parallel: " + time1);
        System.out.println(xi);

        System.out.println("Not Parallel: " + ((end-start)/1000.0));
        System.out.println(xi2);
        //System.out.println(ss.edgeStability(false));
        System.out.println("Directed Edge Stab");
        ss.printEdgeStability(true, new PrintWriter(System.out));
        System.out.println("Undirected Edge Stab");
        ss.printEdgeStability(false, new PrintWriter(System.out));
    }

    public static void tests3(){
        DataSet ds = genTestData();
        double lambda = .05;
        List<Set<Integer>> samps = new ArrayList<>();
        Set<Integer> s1 = new HashSet<Integer>();
        for(int i = 0; i < 100; i++){
            s1.add(i);
        }

        for(int i = 0; i < 1000; i++){
            samps.add(s1);
        }

        //SearchWrappers.MGMWrapper mgm = new SearchWrappers.MGMWrapper(new double[]{lambda, lambda, lambda});
        SearchWrappers.PcStableWrapper mgm = new SearchWrappers.PcStableWrapper(new double[]{.05});
        long start = System.currentTimeMillis();

        //StabilitySearch ss = new StabilitySearch(ds, mgm, N, 500);
        StabilitySearch ss = new StabilitySearch(ds, mgm, samps);
        ss.searchParallel();
        DoubleMatrix2D xi = ss.getThetaMat();
        long end = System.currentTimeMillis();
        double time1 = ((end-start)/1000.0);
        System.out.println("Parallel: " + time1);
        System.out.println(xi);
    }

    public static void listTests(){
        List<Node> variables = new LinkedList<>();
        variables.add(new GraphNode("x1"));
        variables.add(new GraphNode("x2"));
        variables.add(new GraphNode("x3"));

        List<Node> var2 = new LinkedList<>(variables);
        System.out.println("Is x1 equal? " + (variables.get(0)==var2.get(0)));

        //List Set tests
        List<Set<Integer>> test = new ArrayList<Set<Integer>>(3);
        Set<Integer> s1 = new HashSet<Integer>();
        s1.add(1); s1.add(2); s1.add(3);
        test.add(s1);
        Set<Integer> s2 = new HashSet<Integer>();
        s2.add(3); s2.add(1); s2.add(2);
        //test.add(s2);
        Set<Integer> s3 = new HashSet<Integer>();
        s3.add(4); s3.add(5); s3.add(6);
        Set<Integer> s4 = new HashSet<Integer>();
        s4.add(7); s4.add(8); s4.add(9);
        test.add(0,s4);
        System.out.println("Test contains s1?: " + test.contains(s1));
        System.out.println("Test contains s2?: " + test.contains(s2));
        System.out.println("Test contains s3?: " + test.contains(s3));

        //Set<Integer> s1 = new HashSet<Integer>();
        //s1.add(1); s1.add(2); s1.add(3);
        //Set<Integer> c1 = getComplement(s1,7);
        //System.out.println(c1);

        for(Set t : test){
            System.out.println(t);
        }
    }

    //some tests...
    public static void main(String[] args){
        listTests();
        tests3();
    }
}

