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

package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.IndTestMultinomialLogisticRegression;
import edu.cmu.tetrad.search.PcStable;
import edu.cmu.tetrad.search.SearchGraphUtils;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;

/**
 * Created by ajsedgewick on 9/10/15.
 */
public class ExploreIndepTests {
    public static void main(String[] args) {
        testWaldVsMixNew();
        //testThreeVar();
    }

    public static void testThreeVar(){

        try {
            //String path = ExampleMixedSearch.class.getResource("test_data").getPath();
            String path = System.getProperty("user.dir") + "/tetrad-lib/src/main/java/edu/pitt/csb/mgm/test_data";
            System.out.println(path);
            //Graph trueGraph = SearchGraphUtils.patternFromDag(GraphUtils.loadGraphTxt(new File(path, "DAG_0_graph.txt")));
            DataSet dsReal = MixedUtils.loadDataSet(path, "ThreeVarReal.txt");

            IndTestMultinomialAJ indMix = new IndTestMultinomialAJ(dsReal, .05);
            IndTestMultinomialLogisticRegressionWald indWalLin = new IndTestMultinomialLogisticRegressionWald(dsReal, .05, false);
            IndTestMixedMultipleTTest indmtt = new IndTestMixedMultipleTTest(dsReal, .05);

            indMix.setVerbose(true);
            System.out.println(indMix.isIndependent(dsReal.getVariable("X"), dsReal.getVariable("Y"), dsReal.getVariable("Z")));
            System.out.println(indMix.getPValue());

            indmtt.setVerbose(true);
            System.out.println(indmtt.isIndependent(dsReal.getVariable("X"), dsReal.getVariable("Y"), dsReal.getVariable("Z")));
            System.out.println(indmtt.getPValue());
            System.out.println(Arrays.toString(indmtt.dependencePvals(dsReal.getVariable("X"), dsReal.getVariable("Y"), Collections.singletonList(dsReal.getVariable("Z")))));

            PcStableAJ pc = new PcStableAJ(new IndTestMixedMultipleTTest(dsReal, .05));
            pc.setVerbose(true);

            System.out.println(pc.search());

            //IndTestMixedMultipleTTest indWalLog = new IndTestMixedMultipleTTest(ds, .05);

            /*PcStable s1 = new PcStable(indMix);
            PcStable s2 = new PcStable(indWalLin);
            PcStableAJ s3 = new PcStableAJ(indmtt);

            long time = System.currentTimeMillis();
            Graph g1 = SearchGraphUtils.patternFromDag(s1.search());
            System.out.println("Mix Time " + ((System.currentTimeMillis() - time)/1000.0));

            time = System.currentTimeMillis();
            Graph g2 = SearchGraphUtils.patternFromDag(s2.search());
            System.out.println("Wald lin Time " + ((System.currentTimeMillis() - time)/1000.0));

            time = System.currentTimeMillis();
            Graph g3 = SearchGraphUtils.patternFromDag(s3.search());
            System.out.println("Wald log Time " + ((System.currentTimeMillis() - time)/1000.0));

//            System.out.println(g);
//            System.out.println("IndMix: " + s1.search());
//            System.out.println("IndWalLin: " + s2.search());
//            System.out.println("IndWalLog: " + s3.search());

            System.out.println(MixedUtils.EdgeStatHeader);
            System.out.println(MixedUtils.stringFrom2dArray(MixedUtils.allEdgeStats( trueGraph, g1)));
            System.out.println(MixedUtils.stringFrom2dArray(MixedUtils.allEdgeStats( trueGraph, g2)));
            System.out.println(MixedUtils.stringFrom2dArray(MixedUtils.allEdgeStats( trueGraph, g3)));*/
        } catch (Throwable t){
            t.printStackTrace();
        }
    }

    public static void testWaldVsMixNew(){

        try {
            //String path = ExampleMixedSearch.class.getResource("test_data").getPath();
            String path = System.getProperty("user.dir") + "/tetrad-lib/src/main/java/edu/pitt/csb/mgm/test_data";
            System.out.println(path);
            Graph trueGraph = SearchGraphUtils.patternFromDag(GraphUtils.loadGraphTxt(new File(path, "DAG_0_graph.txt")));
            DataSet ds = MixedUtils.loadDataSet(path, "DAG_0_data.txt");

            //IndTestMultinomialAJ indMix = new IndTestMultinomialAJ(ds, .05);
            IndTestFisherZ indz = new IndTestFisherZ(MixedUtils.makeContinuousData(ds), .05);
            //IndTestMultinomialLogisticRegressionWald indWalLin = new IndTestMultinomialLogisticRegressionWald(ds, .05, true);
            IndTestMultinomialLogisticRegression indMLR = new IndTestMultinomialLogisticRegression(ds, .05);
            IndTestMixedMultipleTTest indmtt1 = new IndTestMixedMultipleTTest(MixedUtils.makeContinuousData(ds), .05);
            IndTestMixedMultipleTTest indmtt2 = new IndTestMixedMultipleTTest(MixedUtils.makeContinuousData(ds), .05);
            //IndTestMixedMultipleTTest indWalLog = new IndTestMixedMultipleTTest(ds, .05);

            PcStable s1 = new PcStable(indz);
            PcStable s2 = new PcStable(indmtt1);
            PcStableAJ s3 = new PcStableAJ(indmtt2);
            //s3.setVerbose(true);

            long time = System.currentTimeMillis();
            Graph g1 = SearchGraphUtils.patternFromDag(s1.search());
            System.out.println("Mix Time " + ((System.currentTimeMillis() - time)/1000.0));

            time = System.currentTimeMillis();
            Graph g2 = SearchGraphUtils.patternFromDag(s2.search());
            System.out.println("Wald lin Time " + ((System.currentTimeMillis() - time)/1000.0));

            time = System.currentTimeMillis();
            Graph g3 = SearchGraphUtils.patternFromDag(s3.search());
            System.out.println("Wald log Time " + ((System.currentTimeMillis() - time)/1000.0));

//            System.out.println(g);
//            System.out.println("IndMix: " + s1.search());
//            System.out.println("IndWalLin: " + s2.search());
//            System.out.println("IndWalLog: " + s3.search());

            System.out.println(MixedUtils.EdgeStatHeader);
            System.out.println(MixedUtils.stringFrom2dArray(MixedUtils.allEdgeStats( trueGraph, g1)));
            System.out.println(MixedUtils.stringFrom2dArray(MixedUtils.allEdgeStats( trueGraph, g2)));
            System.out.println(MixedUtils.stringFrom2dArray(MixedUtils.allEdgeStats( trueGraph, g3)));
        } catch (Throwable t){
            t.printStackTrace();
        }
    }

    public static void testWaldVsMix(){
//        Graph g = new EdgeListGraph();
//        g.addNode(new ContinuousVariable("X1"));
//        g.addNode(new ContinuousVariable("X2"));
//        g.addNode(new DiscreteVariable("X3", 4));
//        g.addNode(new DiscreteVariable("X4", 4));
//        g.addNode(new ContinuousVariable("X5"));
//
//        g.addDirectedEdge(g.getNode("X1"), g.getNode("X2"));
//        g.addDirectedEdge(g.getNode("X2"), g.getNode("X3"));
//        g.addDirectedEdge(g.getNode("X3"), g.getNode("X4"));
//        g.addDirectedEdge(g.getNode("X4"), g.getNode("X5"));
//
//        GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(g, "Split(-1.5,-.5,.5,1.5)");
////        System.out.println(pm);
//
//        GeneralizedSemIm im = MixedUtils.GaussianCategoricalIm(pm);
////        System.out.println(im);
//
//        int samps = 200;
//        DataSet ds = im.simulateDataAvoidInfinity(samps, false);
//        ds = MixedUtils.makeMixedData(ds, MixedUtils.getNodeDists(g));
//        //System.out.println(ds);
//        System.out.println(ds.isMixed());
        try {
            //String path = ExampleMixedSearch.class.getResource("test_data").getPath();
            String path = System.getProperty("user.dir") + "/tetrad-lib/src/main/java/edu/pitt/csb/mgm/test_data";
            System.out.println(path);
            Graph trueGraph = SearchGraphUtils.patternFromDag(GraphUtils.loadGraphTxt(new File(path, "DAG_0_graph.txt")));
            DataSet ds = MixedUtils.loadDataSet(path, "DAG_0_data.txt");

            IndTestMultinomialLogisticRegression indMix = new IndTestMultinomialLogisticRegression(ds, .05);
            IndTestMultinomialLogisticRegressionWald indWalLin = new IndTestMultinomialLogisticRegressionWald(ds, .05, true);
            IndTestMultinomialLogisticRegressionWald indWalLog = new IndTestMultinomialLogisticRegressionWald(ds, .05, false);

            PcStable s1 = new PcStable(indMix);
            PcStable s2 = new PcStable(indWalLin);
            PcStable s3 = new PcStable(indWalLog);

            long time = System.currentTimeMillis();
            Graph g1 = SearchGraphUtils.patternFromDag(s1.search());
            System.out.println("Mix Time " + ((System.currentTimeMillis() - time)/1000.0));

             time = System.currentTimeMillis();
            Graph g2 = SearchGraphUtils.patternFromDag(s2.search());
            System.out.println("Wald lin Time " + ((System.currentTimeMillis() - time)/1000.0));

             time = System.currentTimeMillis();
            Graph g3 = SearchGraphUtils.patternFromDag(s3.search());
            System.out.println("Wald log Time " + ((System.currentTimeMillis() - time)/1000.0));

//            System.out.println(g);
//            System.out.println("IndMix: " + s1.search());
//            System.out.println("IndWalLin: " + s2.search());
//            System.out.println("IndWalLog: " + s3.search());

            System.out.println(MixedUtils.EdgeStatHeader);
            System.out.println(MixedUtils.stringFrom2dArray(MixedUtils.allEdgeStats( trueGraph, g1)));
            System.out.println(MixedUtils.stringFrom2dArray(MixedUtils.allEdgeStats( trueGraph, g2)));
            System.out.println(MixedUtils.stringFrom2dArray(MixedUtils.allEdgeStats( trueGraph, g3)));
        } catch (Throwable t){
            t.printStackTrace();
        }
    }

}

