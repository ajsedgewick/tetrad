package edu.pitt.csb.stability;


/**
 * Created by ajsedgewick on 6/5/16.
 */
public abstract class ConstraintGraphSearch extends DataGraphSearch{
    private String testName = "lrt";
    public ConstraintGraphSearch (boolean verbose, double... params){
        super(verbose, params);
    }

    public ConstraintGraphSearch(String test, boolean verbose, double... params){
        super(verbose, params);
        this.testName = test;
    }

    public void setTestName(String name){testName = name;}
    public String getTestName(){return testName;}
}
