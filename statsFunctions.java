/**
 * Utility functions to compute Student's Ttest 
 * http://en.wikipedia.org/wiki/Student's_t-test
 * Apache 2.0 Licensed
 * Author: bzewdu
 */
 
 
public class statsFunctions {

    /**
     * logGamma function
     *
     * @param xx
     * @return
     */
    public static double logGamma(double xx) {
        double stp = 2.5066282746500002D;
        double cof[] = new double[6];
        cof[0] = 76.180091730000001D;
        cof[1] = -86.505320330000004D;
        cof[2] = 24.014098220000001D;
        cof[3] = -1.231739516D;
        cof[4] = 0.00120858003D;
        cof[5] = -5.3638199999999999E-006D;
        double x = xx - 1.0D;
        double tmp = x + 5.5D;
        tmp = (x + 0.5D) * Math.log(tmp) - tmp;
        double ser = 1.0D;
        for (int j = 0; j < 6; j++) {
            x++;
            ser += cof[j] / x;
        }

        return tmp + Math.log(stp * ser);
    }

    /**
     * Gamma function
     *
     * @param x
     * @return
     */
    public static double gamma(double x) {
        double f = 1E+100D;
        double g = 1.0D;
        if (x > 0.0D) {
            for (; x < 3D; x++)
                g *= x;

            f = (1.0D - (2D / (7D * Math.pow(x, 2D))) * (1.0D - 2D / (3D * Math.pow(x, 2D)))) / (30D * Math.pow(x, 2D));
            f = (1.0D - f) / (12D * x) + x * (Math.log(x) - 1.0D);
            f = (Math.exp(f) / g) * Math.pow(6.2831853071795862D / x, 0.5D);
        } else {
            f = (1.0D / 0.0D);
        }
        return f;
    }

    /**
     * utility - betacf
     *
     * @param a
     * @param b
     * @param x
     * @return
     */
    static double betacf(double a, double b, double x) {
        int maxIterations = 50;
        int m = 1;
        double eps = 3.0000000000000001E-005D;
        double am = 1.0D;
        double bm = 1.0D;
        double az = 1.0D;
        double qab = a + b;
        double qap = a + 1.0D;
        double qam = a - 1.0D;
        double bz = 1.0D - (qab * x) / qap;
        for (double aold = 0.0D; m < maxIterations && Math.abs(az - aold) >= eps * Math.abs(az); m++) {
            double tem = (double) m + (double) m;
            double d = ((double) m * (b - (double) m) * x) / ((qam + tem) * (a + tem));
            double ap = az + d * am;
            double bp = bz + d * bm;
            d = (-(a + (double) m) * (qab + (double) m) * x) / ((a + tem) * (qap + tem));
            double app = ap + d * az;
            double bpp = bp + d * bz;
            aold = az;
            am = ap / bpp;
            bm = bp / bpp;
            az = app / bpp;
            bz = 1.0D;
        }

        return az;
    }

    /**
     * betai
     *
     * @param a
     * @param b
     * @param x
     * @return
     */
    public static double betai(double a, double b, double x) {
        double bt = 0.0D;
        double beta;
        if (x == 0.0D || x == 1.0D)
            bt = 0.0D;
        else if (x > 0.0D && x < 1.0D)
            bt = (gamma(a + b) * Math.pow(x, a) * Math.pow(1.0D - x, b)) / (gamma(a) * gamma(b));
        if (x < (a + 1.0D) / (a + b + 2D))
            beta = (bt * betacf(a, b, x)) / a;
        else
            beta = 1.0D - (bt * betacf(b, a, 1.0D - x)) / b;
        return beta;
    }


    /**
     * @param v1
     * @param v2
     * @param f
     * @return
     */
    public static double fDist(double v1, double v2, double f) {
        return betai(v1 / 2D, v2 / 2D, v1 / (v1 + v2 * f));
    }

    /**
     * @param v
     * @return
     */
    public static double student_c(double v) {
        return Math.exp(logGamma((v + 1.0D) / 2D)) / (Math.sqrt(3.1415926535897931D * v) * Math.exp(logGamma(v / 2D)));
    }

    /**
     * @param v
     * @param t
     * @return
     */
    public static double student_tDen(double v, double t) {
        return student_c(v) * Math.pow(1.0D + (t * t) / v, -0.5D * (v + 1.0D));
    }

    /**
     * @param v
     * @param t
     * @return
     */
    public static double stDist(double v, double t) {
        double sm = 0.5D;
        double u = 0.0D;
        double sign = 1.0D;
        double stepSize = t / 5000D;
        if (t < 0.0D)
            sign = -1D;
        else if (t == 0.0D || Double.isInfinite(t) || (Double.isNaN(t)))
            return 0.0D;

        for (u = 0.0D; u <= sign * t; u += stepSize)
            sm += stepSize * student_tDen(v, u);

        if (sign < 0.0D)
            sm = 0.5D - sm;
        else
            sm = 1.0D - sm;
        if (sm < 0.0D)
            sm = 0.0D;
        else if (sm > 1.0D)
            sm = 1.0D;
        return sm;
    }

    /**
     * @param desiredPValue
     * @param am
     * @param av
     * @param an
     * @param bm
     * @param bv
     * @param bn
     * @return
     */
    public static double computePValue(double desiredPValue, double am, double av, int an, double bm, double bv, int bn) {
        _randomVar rv1 = new _randomVar();
        rv1.mean = am;
        rv1.n = an;
        _randomVar rv2 = new _randomVar();
        rv2.mean = bm;
        rv2.n = bn;
        double f = 0;
        double fAlpha = 0;
        if (av > bv) {
            f = av / bv;
            fAlpha = statsFunctions.fDist(rv1.n - 1.0D, rv2.n - 1.0D, f);
        } else {
            f = bv / av;
            fAlpha = statsFunctions.fDist(rv2.n - 1.0D, rv1.n - 1.0D, f);
        }
        double df = 0;
        double t = 0;
        if (fAlpha <= 0.0050000000000000001D) {
            double svn1 = av / rv1.n;
            double svn2 = bv / rv2.n;
            df = Math.pow(svn1 + svn2, 2) / (Math.pow(svn1, 2) / (rv1.n + 1.0D) + Math.pow(svn2, 2) / (rv2.n + 1.0D)) - 2D;
            t = Math.abs(rv1.mean - rv2.mean) / Math.sqrt(av / rv1.n + bv / rv2.n);
        } else {
            df = (rv1.n + rv2.n) - 2D;
            double sp = Math.sqrt(((rv1.n - 1.0D) * av + (rv2.n - 1.0D) * bv) / ((rv1.n + rv2.n) - 2D));
            t = (Math.abs(rv1.mean - rv2.mean) * Math.sqrt((rv1.n * rv2.n) / (rv1.n + rv2.n))) / sp;
        }
        return statsFunctions.stDist(df, t);
    }

    /**
     * @param desiredPValue
     * @param am
     * @param av
     * @param an
     * @param bm
     * @param bv
     * @param bn
     * @return
     */
    public static String computeTtest(double desiredPValue, double am, double av, int an, double bm, double bv, int bn) {
        String comment, comment1;
        _randomVar rv1 = new _randomVar();
        rv1.mean = am;
        rv1.n = an;
        _randomVar rv2 = new _randomVar();
        rv2.mean = bm;
        rv2.n = bn;
        double sv1 = av;
        double sv2 = bv;
        double f = 0;
        double fAlpha = 0;
        statsFunctions sf = new statsFunctions();
        if (sv1 > sv2) {
            f = sv1 / sv2;
            fAlpha = statsFunctions.fDist(rv1.n - 1.0D, rv2.n - 1.0D, f);
        } else {
            f = sv2 / sv1;
            fAlpha = statsFunctions.fDist(rv2.n - 1.0D, rv1.n - 1.0D, f);
        }

        double df = 0;
        double t = 0;
        String newln = "<br>";
        if (fAlpha <= 0.0050000000000000001D) {
            comment = "An F test on the sample variances indicates that they are " + newln + "probably not from the same population (the variances " + newln + "can't be pooled), at an alpha level of " + fAlpha + "." + newln + "Thus, the t-test was set up for samples with unequal varainces. " + newln + "(The degrees of freedom were adjusted.)";
            double svn1 = sv1 / rv1.n;
            double svn2 = sv2 / rv2.n;
            df = Math.pow(svn1 + svn2, 2) / (Math.pow(svn1, 2) / (rv1.n + 1.0D) + Math.pow(svn2, 2) / (rv2.n + 1.0D)) - 2D;
            t = Math.abs(rv1.mean - rv2.mean) / Math.sqrt(sv1 / rv1.n + sv2 / rv2.n);
        } else {
            comment = "An F test on the sample variances indicates that they could be " + newln + "from the same population, (alpha level of 0.005)." + newln + "Accordingly, the t-test was set up for samples with equal population variance.";
            df = (rv1.n + rv2.n) - 2D;
            double sp = Math.sqrt(((rv1.n - 1.0D) * sv1 + (rv2.n - 1.0D) * sv2) / ((rv1.n + rv2.n) - 2D));
            t = (Math.abs(rv1.mean - rv2.mean) * Math.sqrt((rv1.n * rv2.n) / (rv1.n + rv2.n))) / sp;
        }
        double pVal = statsFunctions.stDist(df, t);
        String pValComment = "" + pVal;
        if (pVal <= 0.01D) {
            comment1 = "This probability indicates that there is a difference in sample means.";
            if (pVal <= 0.0001D)
                pValComment = "< 0.0001";
        } else if (pVal <= 0.050000000000000003D)
            comment1 = "This probability indicates that there may be a difference in sample means.";
        else
            comment1 = "There is not a significant difference in the sample means. " + newln +
                    "A difference could not be detected due to large variability, small sample size, " + newln +
                    "or both. Of course, it's possible that the samples really are from the same population!";
        if (rv1.n == 0.0D || rv2.n == 0.0D) {
            comment1 = "There is a problem with the data. Valid delimiters are space, " + newln + "comma, tab and newline.";
            comment = "";
        }
        sv1 = Math.sqrt(sv1);
        sv2 = Math.sqrt(sv2);
        char tab = '\t';
        String cs =
                "Mean of baseline data set                   :" + tab + rv1.mean + newln +
                        "Standard deviation of baseline data set     :" + tab + sv1 + newln +
                        "Number of observations in the baseline set  :" + tab + rv1.n + newln +
                        "Mean of build data set                  :" + tab + rv2.mean + newln +
                        "Standard deviation of build data set    :" + tab + sv2 + newln +
                        "Number of observations in the build set :" + tab + rv2.n + newln + newln +
                        "Degrees of freedom   :" + tab + df + newln +
                        "t Value (one-tailed) :" + tab + t + newln +
                        "P(x>t)               :" + tab + pValComment + newln;
        return cs;
    }


}
