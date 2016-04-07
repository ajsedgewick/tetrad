package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.regression.LogisticRegression;
import edu.cmu.tetrad.regression.LogisticRegressionResult;
import edu.cmu.tetrad.search.IndTestChiSquare;
import edu.cmu.tetrad.search.IndTestGSquare;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.dist.Normal;
import sun.rmi.runtime.Log;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by ajsedgewick on 3/16/16.
 */
public class GeneralTester {

    public static void main(String[] args){
        nullTest();
    }

    private static void nullTest(){
        Node X = new ContinuousVariable("X");
        Node Z = new ContinuousVariable("Z");
        List<Node> vars = new ArrayList<>();
        vars.add(X);
        vars.add(Z);

        int samps = 100;
        int reps = 500;
        Normal n = new Normal(0,1);
        double[][] data = new double[samps][2];
        for (int i = 0; i < samps; i++) {
            data[i][0] = n.nextRandom();
            data[i][1] = n.nextRandom();
        }

        DataSet ds = ColtDataSet.makeData(vars, new TetradMatrix(data));
        IndependenceTest test1 = new IndTestMultinomialAJ(ds, .05);
        IndependenceTest test2 = new IndTestMixedLrt(ds, .05);

        System.out.println("cont TT: " + test1.isIndependent(X, Z) + "\t" + test1.getPValue());
        System.out.println("cont LRT: " + test2.isIndependent(X, Z) + "\t" + test2.getPValue());

    }

    private static void nullTests(){
        Node X = new DiscreteVariable("X",2);
        Node Z = new DiscreteVariable("Z",2);
        List<Node> vars = new ArrayList<>();
        vars.add(X);
        vars.add(Z);

        double zprob = .5;
        double xprob = .5;

        int samps = 100;
        int reps = 500;

        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        int count4 = 0;
        PrintWriter pw = null;
        try {
            pw = new PrintWriter("/Users/ajsedgewick/Desktop/IndTestNull55.txt");
        } catch(FileNotFoundException e) {e.printStackTrace();}

        pw.println("MLL" + "\t" + "TT" + "\t" + "ChiSq" + "\t" + "GSq" + "\t0,0\t1,0\t0,1\t1,1");
        for(int r = 0; r < reps; r++) {
            double[][] data = new double[samps][2];
            for (int i = 0; i < samps; i++) {
                double draw = Math.random();
                if (draw > xprob)
                    data[i][0] = 1.0;
                else
                    data[i][0] = 0.0;

                draw = Math.random();
                if (draw > zprob)
                    data[i][1] = 1.0;
                else
                    data[i][1] = 0.0;
            }


            DataSet ds = ColtDataSet.makeData(vars, new TetradMatrix(data));

            CellTable ct = new CellTable(new int[]{2,2});

            ct.addToTable(ds, new int[]{0,1});

            //System.out.println(ct.calcMargin(new int[]{1,0}));
            //System.out.println(ct);


            //System.out.println(ds.isDiscrete());
            //System.out.println(ds);

            LogisticRegression lr = new LogisticRegression(ds);
            LogisticRegression.Result res = lr.regress((DiscreteVariable) X, Collections.singletonList(Z));
            //System.out.println(res.getLogLikelihood());
            //System.out.println(Arrays.toString(res.getCoefs()));
            //System.out.println(Arrays.toString(res.getStdErrs()));

            IndependenceTest test1 = new IndTestMultinomialAJ(ds, .05);
            //System.out.println("MLL: " + test1.isIndependent(X, Z) + "\t" + test1.getPValue());
            if(!test1.isIndependent(X, Z)) {count1++;}

            IndependenceTest test2 = new IndTestMixedMultipleTTest(ds, .05);
            if(!test2.isIndependent(X, Z)) {count2++;}
            //System.out.println("TT: " + test2.isIndependent(X, Z) + "\t" + test2.getPValue());

            IndependenceTest test3 = new IndTestChiSquare(ds, .05);
            if(!test3.isIndependent(X, Z)) {count3++;}
            //System.out.println("ChiSq: " + test3.isIndependent(X, Z) + "\t" + test3.getPValue());

            IndependenceTest test4 = new IndTestGSquare(ds, .05);
            if(!test4.isIndependent(X, Z)) {count4++;}
            //System.out.println("GSq: " + test4.isIndependent(X, Z) + "\t" + test4.getPValue());
            pw.print(test1.getPValue() + "\t" + test2.getPValue() + "\t" + test3.getPValue() + "\t" + test4.getPValue() + "\t");
            pw.println(ct.calcMargin(new int[]{0,0}) + "\t" + ct.calcMargin(new int[]{1,0}) + "\t" + ct.calcMargin(new int[]{0,1}) + "\t" + ct.calcMargin(new int[]{1,1}));
        }

        //System.out.println("Counts: " + count1 + "\t" + count2 + "\t" + count3 + "\t" + count4);
        pw.flush();
        pw.close();
    }

    private static void runTests3(){
        Graph g = GraphConverter.convert("Z-->Y,Z-->X,Y-->X");
        //Graph g = GraphConverter.convert("X1-->X2");
        //simple graph pm im gen example

        HashMap<String, Integer> nd = new HashMap<>();
        nd.put("X", 0);
        nd.put("Y", 3);
        nd.put("Z", 0);


        g = MixedUtils.makeMixedGraph(g, nd);

        //GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(g, "Split(-1.5,-1,1,1.5)");
        GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(g, "1");
        System.out.println(pm);

        GeneralizedSemIm im = MixedUtils.GaussianCategoricalIm(pm, false);
        System.out.println(im);

        int samps = 100;
        DataSet ds = im.simulateDataAvoidInfinity(samps, false);

        IndependenceTest ind = new IndTestMultinomialAJ(ds, .05);
        System.out.print("Ind(X,Y) | Z: " + ind.isIndependent(ds.getVariable("X"), ds.getVariable("Y"), ds.getVariable("Z")));
        System.out.println(" p: " + Math.log(ind.getPValue()));

        IndTestMixedMultipleTTest ind2 = new IndTestMixedMultipleTTest(ds, .05);
        System.out.print("Ind(X,Y) | Z: " + ind2.isIndependent(ds.getVariable("X"), ds.getVariable("Y"), ds.getVariable("Z")));
        System.out.println(" p: " + Math.log(ind2.getPValue()));

        try{DataWriter.writeRectangularData(ds, new PrintWriter("XYZ_dat.txt"), '\t');}
        catch(Throwable t){t.printStackTrace();}

        DataSet ds2 = DataUtils.getNonparanormalTransformed(ds);

        ind = new IndTestMultinomialAJ(ds2, .05);
        System.out.print("Ind(X,Y) | Z: " + ind.isIndependent(ds2.getVariable("X"), ds2.getVariable("Y"), ds2.getVariable("Z")));
        System.out.println(" p: " + Math.log(ind.getPValue()));

        ind2 = new IndTestMixedMultipleTTest(ds2, .05);
        System.out.print("Ind(X,Y) | Z: " + ind2.isIndependent(ds2.getVariable("X"), ds2.getVariable("Y"), ds2.getVariable("Z")));
        System.out.println(" p: " + Math.log(ind2.getPValue()));


        /*
        ds = MixedUtils.makeMixedData(ds, nd);
        DataSet dsC = MixedUtils.getContinousData(ds);
        DataSet dsD = MixedUtils.getDiscreteData(ds);
        System.out.println("Continuous: rows - " + dsC.getNumRows() + " columns - " + dsC.getNumColumns());
        System.out.println("Discrete: rows - " + dsD.getNumRows() + " columns - " + dsD.getNumColumns());
        //System.out.println(ds);

        double lambda = 0.00000;
        MGM model = new MGM(ds.subsetColumns(new int[]{0}), new double[]{lambda, lambda, lambda});
        //MGM model = new MGM(ds, new double[]{lambda, lambda, lambda});

        System.out.println("Init nll: " + model.smoothValue(model.params.toMatrix1D()));
        System.out.println("Init reg term: " + model.nonSmoothValue(model.params.toMatrix1D()));

        model.learn(1e-8,1000);

        System.out.println("Learned nll: " + model.smoothValue(model.params.toMatrix1D()));
        System.out.println("Learned reg term: " + model.nonSmoothValue(model.params.toMatrix1D()));

        System.out.println("params:\n" + model.params);
        System.out.println("adjMat:\n" + model.adjMatFromMGM());

        DataSet ds2 = im.simulateDataAvoidInfinity(samps, false);
        ds2 = MixedUtils.makeMixedData(ds2, nd);

        //for(Node node : ds.getVariables()){
        Node node = ds.getVariable(0);
        System.out.println("Node: " + node.getName() + " CondNLL: " + model.nodeCondNll(ds.getVariables().indexOf(node)));
        //}

        System.out.println("Var: " + StatUtils.variance(ds.getDoubleData().getColumn(0).toArray()));
        System.out.println("Res: " + .5*ds.getNumRows()* Math.log(StatUtils.variance(ds.getDoubleData().getColumn(0).toArray())/(2*Math.PI)));

        //System.out.println("Test likelihood:\n" + model.smoothValue(model.params.toMatrix1D(), ds2.subsetColumns(new int[]{1})));
        //System.out.println("Test likelihood:\n" + model.smoothValue(model.params.toMatrix1D(), ds2));
        */
    }
}
