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
import org.apache.commons.math3.util.CombinatoricsUtils;

import java.io.File;
import java.io.PrintStream;
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
    final private int[][] samps;
    final private DoubleMatrix2D thetaMat;
    final private HashMap<Edge, Double> edgeStab = new HashMap<>();
    final private DataGraphSearch gs;
    private int numVars;
    private int numSamps;
    private int numSubs;
    private int subSize;

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

        this.samps = subSampleNoReplacement(numSamps, numSubs, subSize);
        this.thetaMat = DoubleFactory2D.dense.make(numVars, numVars, 0.0);
        //this(data, gs, 20, subSize);
    }

    public StabilitySearch(DataSet data, DataGraphSearch gs, int N, int b){
        this.data = data;
        this.gs = gs;
        this.numSamps = data.getNumRows();
        this.numSubs = N;
        this.numVars = data.getNumColumns();
        this.subSize = b;

        this.samps = subSampleNoReplacement(numSamps, N, b);
        this.thetaMat = DoubleFactory2D.dense.make(numVars, numVars, 0.0);
    }

    //returns an adjacency matrix containing the edgewise instability as defined in Liu et al
    public DoubleMatrix2D searchSerial(){

        for(int s = 0; s < numSubs; s++){
            DataSet dataSubSamp = data.subsetRows(samps[s]);
            Graph g = gs.search(dataSubSamp);

            DoubleMatrix2D curAdj = MixedUtils.skeletonToMatrix(g);
            thetaMat.assign(curAdj, Functions.plus);
        }
        thetaMat.assign(Functions.mult(1.0 / numSubs));
        //thetaMat.assign(thetaMat.copy().assign(Functions.minus(1.0)), Functions.mult).assign(Functions.mult(-2.0));
        return thetaMat;
    }

    //returns an adjacency matrix containing the edgewise instability as defined in Liu et al
    public  DoubleMatrix2D searchParallel(){

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
            private void addToMat(DoubleMatrix2D matSum, DoubleMatrix2D curMat){
                matSum.assign(curMat, Functions.plus);
                if(!matSum.equals(Algebra.DEFAULT.transpose(matSum))) {
                    System.out.println("NOT SYMMETRIC!!!:\n" + matSum + "\ncurmat:\n" + curMat);
                    //throw new IllegalStateException("not symmetric");
                }
            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        DataSet dataSubSamp = data.subsetRows(samps[s]).copy();
                        DataGraphSearch curGs = gs.copy();
                        Graph g = curGs.search(dataSubSamp);

                        synchronized (edgeStab){
                            addEdges(g);
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

        final int chunk = 1;

        pool.invoke(new StabilityAction(chunk, 0, numSubs));

        thetaMat.assign(Functions.mult(1.0 / numSubs));

        //do this elsewhere
        //thetaMat.assign(thetaMat.copy().assign(Functions.minus(1.0)), Functions.mult).assign(Functions.mult(-2.0));
        return thetaMat;
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

    /**
     * print list of edge stabilities to tab-delimited format, round to 5 decimal places. currently does not handle
     * circle endpoints...
     * @param directed
     * @param ps
     */
    public void printEdgeStability(boolean directed, PrintStream ps){
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

            ps.println(e.getNode1() + dlm + e.getNode2() + dlm + eType + dlm + df.format(d));
        }
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
    public static int[][] subSampleNoReplacement(int sampSize, int numSub, int subSize){

        if (subSize < 1) {
            throw new IllegalArgumentException("Sample size must be > 0.");
        }

        List<Integer> indices = new ArrayList<Integer>(sampSize);
        for (int i = 0; i < sampSize; i++) {
            indices.add(i);
        }

        int[][] sampMat = new int[numSub][subSize];

        for(int i = 0; i < numSub; i++) {
            Collections.shuffle(indices);
            int[] curSamp;
            SAMP:
            while(true){
                curSamp = subSampleIndices(sampSize, subSize);
                for(int j = 0; j < i; j++){
                    if(Arrays.equals(curSamp, sampMat[j])){
                        continue SAMP;
                    }
                }
                break;
            }
            sampMat[i] = curSamp;
        }

        return sampMat;
    }

    public static int[] subSampleIndices(int N, int subSize){
        List<Integer> indices = new ArrayList<Integer>(N);
        for (int i = 0; i < N; i++) {
            indices.add(i);
        }

        Collections.shuffle(indices);
        int[] samp = new int[subSize];
        for(int i = 0; i < subSize; i++){
            samp[i] = indices.get(i);
        }
        return samp;
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
        DoubleMatrix2D xi = ss.searchSerial();
        long end = System.currentTimeMillis();
        System.out.println("Not parallel: " + ((end-start)/1000.0));

        start = System.currentTimeMillis();
        ss = new StabilitySearch(ds, mgm, 20, 200);
        DoubleMatrix2D xi2 = ss.searchParallel();
        end = System.currentTimeMillis();
        System.out.println("Parallel: " + ((end-start)/1000.0));

        System.out.println(xi);
        System.out.println(xi2);
    }

    private static void tests2(){
        Graph g = GraphConverter.convert("X1-->X2,X3-->X2,X4-->X5");
        //Graph g = GraphConverter.convert("X1-->X2");
        //simple graph pm im gen example

        HashMap<String, Integer> nd = new HashMap<>();
        nd.put("X1", 0);
        nd.put("X2", 3);
        nd.put("X3", 3);
        nd.put("X4", 3);
        nd.put("X5", 0);

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

        double lambda = .05;
        //SearchWrappers.MGMWrapper mgm = new SearchWrappers.MGMWrapper(new double[]{lambda, lambda, lambda});
        SearchWrappers.PcStableWrapper mgm = new SearchWrappers.PcStableWrapper(new double[]{.05});
        long start = System.currentTimeMillis();
        int N = 20;
        StabilitySearch ss = new StabilitySearch(ds, mgm, N, 500);
        DoubleMatrix2D xi = ss.searchSerial();
        long end = System.currentTimeMillis();
        double time1 = ((end-start)/1000.0);

        start = System.currentTimeMillis();
        ss = new StabilitySearch(ds, mgm, N, 500);
        DoubleMatrix2D xi2 = ss.searchParallel();
        end = System.currentTimeMillis();

        System.out.println("Not parallel: " + time1);
        System.out.println(xi);

        System.out.println("Parallel: " + ((end-start)/1000.0));
        System.out.println(xi2);
        //System.out.println(ss.edgeStability(false));
        System.out.println("Directed Edge Stab");
        ss.printEdgeStability(true, System.out);
        System.out.println("Undirected Edge Stab");
        ss.printEdgeStability(false, System.out);
    }



    //some tests...
    public static void main(String[] args){
        tests2();
    }
}

