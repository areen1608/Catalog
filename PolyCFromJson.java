import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Unique solution:
 * - Custom base decoder (no BigInteger(value, radix))
 * - Barycentric Lagrange at x=0 (exact rational), c = f(0)
 * - Arbitrary precision via BigInteger + minimal BigRational
 *
 * Usage (with Jackson jars on classpath):
 * javac -cp
 * .:jackson-databind-2.17.2.jar:jackson-core-2.17.2.jar:jackson-annotations-2.17.2.jar
 * PolyCFromJson.java
 * java -cp
 * .:jackson-databind-2.17.2.jar:jackson-core-2.17.2.jar:jackson-annotations-2.17.2.jar
 * PolyCFromJson testcase1.json
 */
public class PolyCFromJson {

    // immutable pair
    static final class XY {
        final BigInteger x;
        final BigInteger y;

        XY(BigInteger x, BigInteger y) {
            this.x = x;
            this.y = y;
        }
    }

    // minimal exact rational
    static final class BigRational {
        final BigInteger num; // can be negative
        final BigInteger den; // > 0

        BigRational(BigInteger n, BigInteger d) {
            if (d.signum() == 0)
                throw new ArithmeticException("denominator = 0");
            if (d.signum() < 0) {
                n = n.negate();
                d = d.negate();
            }
            BigInteger g = n.gcd(d);
            this.num = g.equals(BigInteger.ONE) ? n : n.divide(g);
            this.den = g.equals(BigInteger.ONE) ? d : d.divide(g);
        }

        static BigRational of(BigInteger z) {
            return new BigRational(z, BigInteger.ONE);
        }

        BigRational add(BigRational r) {
            return new BigRational(num.multiply(r.den).add(r.num.multiply(den)), den.multiply(r.den));
        }

        BigRational sub(BigRational r) {
            return new BigRational(num.multiply(r.den).subtract(r.num.multiply(den)), den.multiply(r.den));
        }

        BigRational mul(BigRational r) {
            return new BigRational(num.multiply(r.num), den.multiply(r.den));
        }

        BigRational div(BigRational r) {
            if (r.num.signum() == 0)
                throw new ArithmeticException("divide by zero");
            return new BigRational(num.multiply(r.den), den.multiply(r.num));
        }

        boolean isInteger() {
            return den.equals(BigInteger.ONE);
        }

        @Override
        public String toString() {
            return isInteger() ? num.toString() : (num + "/" + den);
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.err.println("usage: java PolyCFromJson <path/to/testcase.json>");
            System.exit(1);
        }
        String raw = Files.readString(Paths.get(args[0]));
        ObjectMapper om = new ObjectMapper();
        JsonNode root = om.readTree(raw);

        // read k (and n, though we only need k for degree)
        JsonNode keys = must(root, "keys");
        int n = must(keys, "n").asInt();
        int k = must(keys, "k").asInt();
        if (k < 1)
            fail("k must be >= 1");
        if (n < k)
            fail("n must be >= k");

        // collect (x, y) decoding y with our custom parser
        List<XY> points = new ArrayList<>();
        Iterator<String> it = root.fieldNames();
        while (it.hasNext()) {
            String fname = it.next();
            if ("keys".equals(fname))
                continue;
            if (!fname.chars().allMatch(Character::isDigit))
                continue; // numeric keys only

            BigInteger x = new BigInteger(fname);
            JsonNode obj = root.get(fname);

            String baseStr = must(obj, "base").asText();
            String valStr = must(obj, "value").asText();
            int base = parseBase(baseStr);
            BigInteger y = parseInBase(valStr, base);

            points.add(new XY(x, y));
        }

        if (points.size() < k)
            fail("need at least k=" + k + " points in JSON");
        // sort by x, then take first k distinct x
        points.sort(Comparator.comparing(p -> p.x));
        List<XY> chosen = new ArrayList<>();
        for (XY p : points) {
            if (chosen.stream().noneMatch(q -> q.x.equals(p.x))) {
                chosen.add(p);
                if (chosen.size() == k)
                    break;
            }
        }
        if (chosen.size() < k)
            fail("not enough distinct x's to reach k=" + k);

        // compute c = f(0) using barycentric form (exact)
        BigRational c = interpolateAtZeroBarycentric(chosen);

        // output
        System.out.println("k = " + k + " (degree " + (k - 1) + ")");
        System.out.println("points used:");
        for (XY p : chosen)
            System.out.println("  x=" + p.x + ", y=" + p.y);
        System.out.println("\nmethod: barycentric lagrange @ x=0 (exact)");
        System.out.println("secret c = " + c.toString());
    }

    // ---------- math: barycentric lagrange at x = 0 ----------
    // f(0) = [sum_i (w_i * y_i / -x_i)] / [sum_i (w_i / -x_i)] = [sum_i
    // (w_i*y_i/x_i)] / [sum_i (w_i/x_i)]
    static BigRational interpolateAtZeroBarycentric(List<XY> pts) {
        int k = pts.size();
        // precompute barycentric weights w_i = 1 / Π_{j≠i} (x_i - x_j)
        BigRational[] w = new BigRational[k];
        for (int i = 0; i < k; i++) {
            BigInteger denom = BigInteger.ONE;
            for (int j = 0; j < k; j++) {
                if (i == j)
                    continue;
                BigInteger diff = pts.get(i).x.subtract(pts.get(j).x);
                if (diff.signum() == 0)
                    fail("duplicate x encountered");
                denom = denom.multiply(diff);
            }
            w[i] = new BigRational(BigInteger.ONE, denom);
        }

        BigRational num = BigRational.of(BigInteger.ZERO);
        BigRational den = BigRational.of(BigInteger.ZERO);
        for (int i = 0; i < k; i++) {
            XY pi = pts.get(i);
            if (pi.x.signum() == 0)
                fail("x_i = 0 not supported in barycentric at x=0");
            BigRational invXi = new BigRational(BigInteger.ONE, pi.x); // 1/x_i
            BigRational wy = w[i].mul(BigRational.of(pi.y)); // w_i * y_i
            num = num.add(wy.mul(invXi)); // sum w_i*y_i/x_i
            den = den.add(w[i].mul(invXi)); // sum w_i/x_i
        }
        return num.div(den);
    }

    // ---------- JSON helpers ----------
    private static JsonNode must(JsonNode node, String key) {
        if (node == null || !node.has(key))
            fail("missing field: " + key);
        return node.get(key);
    }

    private static void fail(String msg) {
        System.err.println(msg);
        System.exit(2);
    }

    // ---------- unique base parser (no BigInteger(value, radix)) ----------
    private static int parseBase(String s) {
        try {
            int b = Integer.parseInt(s.trim());
            if (b < 2 || b > 36)
                fail("base out of range (2..36): " + s);
            return b;
        } catch (NumberFormatException e) {
            fail("invalid base: " + s);
            return -1; // unreachable
        }
    }

    private static BigInteger parseInBase(String digits, int base) {
        BigInteger B = BigInteger.valueOf(base);
        BigInteger acc = BigInteger.ZERO;
        for (int i = 0; i < digits.length(); i++) {
            int v = valueOfDigit(digits.charAt(i));
            if (v < 0 || v >= base)
                fail("digit '" + digits.charAt(i) + "' not valid in base " + base);
            acc = acc.multiply(B).add(BigInteger.valueOf(v));
        }
        return acc;
    }

    private static int valueOfDigit(char ch) {
        if (ch >= '0' && ch <= '9')
            return ch - '0';
        if (ch >= 'a' && ch <= 'z')
            return 10 + (ch - 'a');
        if (ch >= 'A' && ch <= 'Z')
            return 10 + (ch - 'A');
        fail("invalid digit: " + ch);
        return -1;
    }
}
