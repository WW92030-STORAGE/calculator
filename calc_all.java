package ALLINONE;

import java.util.ArrayList;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;

class info { // meta info (error messages and help)
	static String[] messages = {"SYNTAX ERROR", "SYNTAX ERROR", "SYNTAX ERROR", "SYNTAX ERROR", "SYNTAX ERROR", "SYNTAX ERROR", "SYNTAX ERROR", "SYNTAX ERROR", "SYNTAX ERROR", "bruh"};
	static TreeSet<String> errors = new TreeSet<String>();
	static TreeSet<String> weird = new TreeSet<String>();
	
	static void init() {
		System.out.println("CALCULATOR");
		System.out.println("Type help or functions for more information");
		
		calc_all.ans = 0;
		errors.clear();
		for (String s : messages) errors.add(s.toLowerCase());
		
		weird.add("calc");
		weird.add("calculator");
	}
	
	static void error() {
		int rng = (int)(Math.random() * messages.length);
		System.out.println(messages[rng]);
	}
	
	static void help() {
		System.out.println("Type in an expression and the computer will magically solve it\nUse ans to reference previous results\nNumbers above 10^308 are not supported");
	}
	
	static void asset() {
		System.out.println("SUPPORTED FUNCTIONS (all functions except for rad are in radian mode)");
		System.out.println("sin cos tan sec csc cot rad exp ln log phi inv sigma tau mu");
		System.out.println("Type help <function> for information on <function>");
	}
	
	static void help(String e) {
				// TRIG
				if (e.equals("sin")) System.out.println("sin(x) - returns the sine in radians"); 
				else if (e.equals("cos")) System.out.println("cos(x) - returns the cosine in radians"); 
				else if (e.equals("tan")) System.out.println("tan(x) - returns the tangent in radians");
				else if (e.equals("csc")) System.out.println("csc(x) - returns the cosecant in radians"); 
				else if (e.equals("sec")) System.out.println("sec(x) - returns the secant in radians"); 
				else if (e.equals("cot")) System.out.println("cot(x) - returns the cotangent in radians"); 
				
				// NUMBER THEORY
				else if (e.equals("phi")) System.out.println("phi(n) - returns the number of positive integers below n relatively prime to floor(n)"); 
				else if (e.equals("exp")) System.out.println("exp(b, e) - returns the result b^e"); 
				else if (e.equals("inv")) System.out.println("inv(a, m) - returns the modular inverse of a mod floor(m)\nuse at your own risk with non-integer values of a\nif gcd(a, m) != 1 then the function pulls a cheese move..."); 
				else if (e.equals("log")) System.out.println("log(b, n) - returns the logarithm base b of n"); 
				else if (e.equals("ln")) System.out.println("ln(n) - returns the logarithm base e of n"); 
				else if (e.equals("tau")) System.out.println("tau(n) - returns the number of positive integer divisors of n"); 
				else if (e.equals("rad")) System.out.println("rad(n) - converts n from degrees to radians"); 
				else if (e.equals("sigma")) System.out.println("sigma(n) - returns the sum of the postive integer divisors of n"); 
				else if (e.equals("mu")) System.out.println("mu(n) - returns 0 if has square divisor\nreturns 1 if no squared divisor and even number of prime divisors\nreturns -1 if no squared divisor and odd number of prime divisors");
				else if (e.equals("gcd")) System.out.println("gcd(a, b) - returns the greatest common divisor of a and b\nWARNING - impractical for irrational numbers");
				
				// ARITHMETIC
				else if (e.equals("sqrt")) System.out.println("sqrt(n) - gives the positive root of the equation x^2 = n");
				
				else System.out.println("INVALID NAME");
	}
	
	static void idk() {
		long start = System.nanoTime();
		
		int limit = 4;
		
		for (int i = 0; i < limit;) {
			if (System.nanoTime() > start + 500 * 1000 * 1000) {
				if (i < limit - 1) System.out.print("...");
				start = System.nanoTime();
				i++;
			}
			
		}
		
		System.out.print("why\n");
	}
}

class functions { // functions library
	// UTILITY
	public static int factorial(int x) {
		if (x <= 1) return 1;
		int res = 1;
		for (int i = 2; i <= x; i++) res *= i;
		return res;
	}
	
	public static double exp(double b, int e) {
		if (e == 0) return 1;
		double half = exp(b, e / 2);
		if (e % 2 == 0) return half * half;
		return half * half * b;
	}
	
	public static double exp(double b, int e, double m) {
		if (e == 0) return 1;
		double half = exp(b, e / 2, m);
		if (e % 2 == 0) return (half * half) % m;
		return (half * half * b) % m;
	}
	
	public static int exp(int b, int e, int m) {
		if (e == 0) return 1;
		int half = exp(b, e / 2, m);
		if (e % 2 == 0) return (half * half) % m;
		return (half * half * b) % m;
	}
	
	public static TreeSet<Double> factorize(double n) {
		double x = n;
		TreeSet<Double> factors = new TreeSet<Double>();
		double i = 2;
		while (i <= x) {
			if (x % i == 0) {
				factors.add(i);
				while (x % i == 0) x /= i;
			}
			i++;
		}
		
		return factors;
	}
	
	public static TreeMap<Double, Double> expansion(double n) {
		double x = n;
		TreeMap<Double, Double> factors = new TreeMap<Double, Double>();
		double i = 2;
		while (i <= x) {
			if (x % i == 0) {
				factors.put(i, 0.0);
				while (x % i == 0) {
					x /= i;
					factors.put(i, factors.get(i) + 1);
				}
			}
			i++;
		}
		
		return factors;
	}
	
	public static String printExpansion(double n) { 
		TreeMap<Double, Double> factors = expansion(n);
		
		String s = "";
		for (double i : factors.keySet()) {
			if (s.length() > 0) s = s + " * ";
			s = s + "" + i + "^" + factors.get(i);
		}
		return s;
	}
	
	//FUNCTIONS
	//READ THE DESCRIPTIONS IN THE INFO CLASS
	
	
	public static double phi(double x) {
		TreeSet<Double> n = factorize(x);
		
		for (double i : n) {
			x *= (i - 1);
			x /= i;
		}
		
		return x;
	}
	
	public static double inv(double a, int m) {
		if (gcd(a, m) != 1.0) return 1.0 / a;
		double p = phi(m);
		return exp(a, (int)p - 1, m);
	}
	
	public static double tau(int x) {
		TreeMap<Double, Double> stuff = expansion(x);
		double res = 1;
		for (double i : stuff.keySet()) {
			res *= (stuff.get(i) + 1);
		}
		return res;
	}
	
	public static double sigma(int x) {
		TreeMap<Double, Double> stuff = expansion(x);
		double res = 1;
		for (double i : stuff.keySet()) {
			double top = Math.pow(i, stuff.get(i) + 1) - 1;
			double bottom = i - 1;
			res *= (top / bottom);
		}
		return res;
	}
	
	public static double mu(int x) {
		TreeMap<Double, Double> stuff = expansion(x);
		for (double i : stuff.keySet()) if (stuff.get(i) >= 2) return 0;
		
		return exp(-1, stuff.size());
	}
	
	public static double gcd(double a, double b) {
		double epsilon = Math.pow(10, -10);
		if (a < b) {
			double temp = a;
			a = b;
			b = temp;
		}
		if (Math.abs(b) < epsilon) return a;
		
		return gcd(b, a % b);
	}
}

class evaluate { // quick maths
	static ArrayList<String> errorMessage = new ArrayList<String>(); // actually this is unused
	static boolean isNumber(String x) {
		for (int i = 0; i < x.length(); i++) {
			if (!parser.isDigit(x.charAt(i))) return false;
		}
		return true;
	}
	
	static ArrayList<String> manage(ArrayList<String> expressio) { // clean up 
		ArrayList<String> expression = new ArrayList<String>();
		
		
		
		for (int i = 0; i < expressio.size(); i++) { // multiplication within parentheses e.g. 3(2 + 4) = 3 * (2 + 4)
			String s = expressio.get(i);
			expression.add(s);
			if (i < expressio.size() - 1) {
				String next = expressio.get(i + 1);
				boolean closednum = s.equals(")") && isNumber(next); // (3)3
				boolean twonum = isNumber(s) && isNumber(next); // 3 [*] 3
				boolean paren = s.equals(")") && next.equals("("); // (3)(3)
				boolean numopen = isNumber(s) && next.equals("("); // 3(3) 
				boolean parenfunc = s.equals(")") && parser.isFunction(next); // (3) ln(3) 
				boolean numfunc = isNumber(s) && parser.isFunction(next); // 3ln(3)
				if (closednum || twonum || paren || numopen || parenfunc || numfunc) expression.add("*");
			}
		}
		
		// note that for parenfunc if the function was before the number or another function then it would fall into closednum
		
		if (expression.get(0).equals("-")) { // negative sign at beginning
			expression.remove(0);
			expression.add(0, "*");
			expression.add(0, "-1");
		}
		
		for (int i = 1; i < expression.size(); i++) {
			String now = expression.get(i);
			String prev = expression.get(i - 1);
			if (now.equals("-")) { // manage the negative sign
				if (prev.equals("+")) {
					i--;
					expression.remove(i);
				}
				else if (prev.equals("-")) { // 
					i--;
					expression.remove(i);
					expression.set(i, "+");
				}
				else if (prev.equals("*")) { // multiplication by negative number
					expression.set(i, "-1");
					expression.add(i + 1, "*");
				}
				else if (prev.equals("/") && !parser.isFunction(expression.get(i + 1))) { // wait until the function is resolved and then apply negatives
					double manip = Double.parseDouble(expression.get(i + 1));
					expression.set(i, "-1");
					expression.add(i + 1, "*");
					expression.set(i + 2, Double.toString(1 / manip));
				}
			}
		}
		
		return expression;
	}
	
	static double evaluate(ArrayList<String> expression) {
		String p = Double.toString(Math.PI);
		String ee = Double.toString(Math.E);
		for (int i = 0; i < expression.size(); i++) {
			String s = expression.get(i);
			if (s.equals("pi")) expression.set(i, p);
			else if (s.equals("e")) expression.set(i, ee);
			else if (s.equals("ans")) expression.set(i, Double.toString(calc_all.ans));
		}
		
		if (expression.get(0).equals("-")) {
			expression.set(0, "-1");
			expression.add(1, "*");
		}
		
	//	System.out.println(expression);
		expression = bracket(expression); // FOR DESCRIPTIONS OF THE METHODS SEE THE evalInside METHOD
		
	//	System.out.println(expression);
		
		expression = manage(expression);
		
	//	System.out.println(expression);
		
		expression = func(expression);
		expression = manage(expression);
		expression = exp(expression);
		expression = manage(expression);
		expression = md(expression);
		expression = manage(expression);
		expression = as(expression);
		expression = manage(expression);
	//	System.out.println(expression);
		return Double.parseDouble(expression.get(0));
	}

	static ArrayList<String> evalInside(ArrayList<String> expression) { // READ FOR DESCRIPTIONS
	//	System.out.println(expression);
		expression = bracket(expression); // recursively resolve bracketed groups
		
		expression = manage(expression); // manage simply makes notation readable by the computer		
	//	System.out.println(expression);
		
		expression = func(expression); // resolve functions
		expression = manage(expression);
		expression = exp(expression); // by this point there should only be numbers and raw operations
		expression = manage(expression);
		expression = md(expression); // multiplication stuff
		expression = manage(expression);
		expression = as(expression); // addition stuff
		expression = manage(expression);
	//	System.out.println(expression);
		return expression; // slightly different
	}
	
	static ArrayList<String> bracket(ArrayList<String> expression) { // recursively simpify stuff in the brackets
	//	System.out.println(expression);
		
		int curLevel = 0;
		int openLevel = -1;
		int open = -1;
		for (int i = 0; i < expression.size(); i++) {
			
			if (expression.get(i).equals("(")) { // found open parentheses
				
				curLevel++;
				if (openLevel == -1) {
					openLevel = curLevel;
					open = i;
				}
				
			//	System.out.println("OPEN " + curLevel + " " + i + " " + openLevel);
			}
		//	if (expression.get(i).equals(")")) System.out.println("CLOSE " + curLevel + " " + i);
			if (expression.get(i).equals(")") && curLevel == openLevel) { // first closing parentheses on same level
			//	System.out.println("FOUND MATCH");
				ArrayList<String> inner = new ArrayList<String>();
				for (int j = open + 1; j < i; j++) {
					String s = expression.get(open + 1);
					inner.add(s);
					expression.remove(open + 1);
				}
				
			//	System.out.println("INNER LEVEL " + openLevel + " : " + inner);
				ArrayList<String> res = evalInside(inner);
			//	System.out.println("INNER = " + res);
				for (int k = res.size() - 1; k >= 0; k--) expression.add(open + 1, res.get(k));
				expression.remove(open + res.size() + 1);
				expression.remove(open);
				i = open;
				openLevel = -1;
			}
			if (expression.get(i).equals(")")) {
				curLevel--;
			}
		}
		
	//	System.out.println(expression);
		return expression;
	}
	
	// list of functions
	// sin cos tan csc sec cot
	// phi exp inv log ln sqrt
	
	static ArrayList<String> func(ArrayList<String> expression) {
		for (int i = 0; i < expression.size(); i++) {
			String s = expression.get(i);
			if (parser.isFunction(s)) {
				double next = Double.parseDouble(expression.get(i + 1));
				double twice;
				if (s.equals("sin")) {
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(Math.sin(next)));
				}
				if (s.equals("cos")) {
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(Math.cos(next)));
				}
				if (s.equals("tan")) {
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(Math.tan(next)));
				}
				if (s.equals("csc")) {
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(1 / Math.sin(next)));
				}
				if (s.equals("sec")) {
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(1 / Math.cos(next)));
				}
				if (s.equals("cot")) {
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(1 / Math.tan(next)));
				}
				if (s.equals("rad")) {
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(next * Math.PI / 180.0));
				}
				if (s.equals("exp")) {
					twice = Double.parseDouble(expression.get(i + 3)); // remember comma separates parameters
					expression.remove(i);
					expression.remove(i);
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(Math.pow(next, twice)));
				}
				if (s.equals("log")) {
					twice = Double.parseDouble(expression.get(i + 3)); // remember comma separates parameters
					expression.remove(i);
					expression.remove(i);
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(Math.log(twice) / Math.log(next)));
				}
				if (s.equals("ln")) {
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(Math.log(next)));
				}
				if (s.equals("inv")) {
					twice = Double.parseDouble(expression.get(i + 3)); // remember comma separates parameters
					expression.remove(i);
					expression.remove(i);
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(functions.inv(next, (int)(twice))));
				}
				if (s.equals("phi")) { // round to int first
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(functions.phi((int)(next))));
				}
				if (s.equals("tau")) { // round to int first
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(functions.tau((int)(next))));
				}
				if (s.equals("sigma")) { // round to int first
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(functions.sigma((int)(next))));
				}
				if (s.equals("mu")) { // round to int first
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(functions.mu((int)(next))));
				}
				if (s.equals("gcd")) {
					twice = Double.parseDouble(expression.get(i + 3)); // remember comma separates parameters
					expression.remove(i);
					expression.remove(i);
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(functions.gcd(next, twice)));
				}
				
				if (s.equals("sqrt")) {
					expression.remove(i);
					expression.remove(i);
					expression.add(i, Double.toString(Math.sqrt(next)));
				}
			}
			
		}
		
		return expression;
	}
	
	
	
	
	
	static ArrayList<String> exp(ArrayList<String> expression) { // ^ 
		for (int i = 0; i < expression.size(); i++) {
			String s = expression.get(i);
			if (s.equals("^")) {
				double first = Double.parseDouble(expression.get(i - 1));
				double second = Double.parseDouble(expression.get(i + 1));
				for (int j = 0; j < 3; j++) expression.remove(i - 1);
				expression.add(i - 1, Double.toString(Math.pow(first, second)));
				i--;
			}
		}
		return expression;
	}
	
	static ArrayList<String> md(ArrayList<String> expression) { // */%
		for (int i = 0; i < expression.size(); i++) {
			String s = expression.get(i);
			if (s.equals("*")) {
				double first = Double.parseDouble(expression.get(i - 1));
				double second = Double.parseDouble(expression.get(i + 1));
				for (int j = 0; j < 3; j++) expression.remove(i - 1);
				expression.add(i - 1, Double.toString(first * second));
				i--;
			}
			
			if (s.equals("/")) {
				double first = Double.parseDouble(expression.get(i - 1));
				double second = Double.parseDouble(expression.get(i + 1));
				for (int j = 0; j < 3; j++) expression.remove(i - 1);
				expression.add(i - 1, Double.toString(first / second));
				i--;
			}
			
			if (s.equals("%")) {
				double first = Double.parseDouble(expression.get(i - 1));
				double second = Double.parseDouble(expression.get(i + 1));
				for (int j = 0; j < 3; j++) expression.remove(i - 1);
				expression.add(i - 1, Double.toString(first % second));
				i--;
			}
		}
		return expression;
	}
	
	static ArrayList<String> as(ArrayList<String> expression) { // +-
		for (int i = 0; i < expression.size(); i++) {
			String s = expression.get(i);
			if (s.equals("+")) {
				double first = Double.parseDouble(expression.get(i - 1));
				double second = Double.parseDouble(expression.get(i + 1));
				for (int j = 0; j < 3; j++) expression.remove(i - 1);
				expression.add(i - 1, Double.toString(first + second));
				i--;
			}
			else if (s.equals("-")) {
				double first = Double.parseDouble(expression.get(i - 1));
				double second = Double.parseDouble(expression.get(i + 1));
				for (int j = 0; j < 3; j++) expression.remove(i - 1);
				expression.add(i - 1, Double.toString(first - second));
				i--;
			}
		}
		return expression;
	}
	
}

class parser { // PARSE ONLY
	
	static boolean isFunction(String e) {
		// TRIG
		if (e.equals("sin")) return true; 
		if (e.equals("cos")) return true;
		if (e.equals("tan")) return true;
		if (e.equals("csc")) return true;
		if (e.equals("sec")) return true;
		if (e.equals("cot")) return true;
		
		// NUMBER THEORY
		if (e.equals("phi")) return true;
		if (e.equals("exp")) return true;
		if (e.equals("inv")) return true;
		if (e.equals("log")) return true;
		if (e.equals("ln")) return true;
		if (e.equals("tau")) return true;
		if (e.equals("rad")) return true;
		if (e.equals("sigma")) return true;
		if (e.equals("mu")) return true;
		if (e.equals("gcd")) return true;
		
		// ARITHMETIC
		if (e.equals("sqrt")) return true;
		
		return false;
	}
	
	static boolean isConstant(String e) { 
		if (e.equals("pi")) return true;
		if (e.equals("e")) return true;
		if (e.equals("ans")) return true;
		return false;
	}
	
	static boolean isOperation(char c) {
		if (c == '+') return true;
		if (c == '-') return true;
		if (c == '*') return true;
		if (c == '/') return true;
		if (c == '%') return true;
		if (c == '^') return true;
		if (c == '(') return true;
		if (c == ')') return true;
		return false;
	}
	
	static boolean isDigit(char c) {
		if (c >= '0' && c <= '9') return true;
		return c == '.';
	}
	
	static boolean isLetter(char c) {
		if (c >= 'A' && c <= 'Z') return true;
		return c >= 'a' && c <= 'z';
	}
	
	static ArrayList<String> parse(String e) {
		for (int i = 0; i < 4; i++) e = e + " ";
		ArrayList<String> parsed = new ArrayList<String>();
		String running = "";
		for (int i = 0; i < e.length(); i++) {
			char c = e.charAt(i);
		//	System.out.println(">> " + c + " <<");
			if (c == ',') {
			//	System.out.println("COMMA");
				parsed.add(running);
				parsed.add(",");
				running = "";
			}
			else if (c == ' ') {
			//	System.out.println("WHITESPACE");
				parsed.add(running);
				running = "";
			}
			else if (isOperation(c)) {
			//	System.out.println("OPERATION");
				parsed.add(running);
				running = "";
				parsed.add("" + c);
			}
			else if (isDigit(c)) {
			//	System.out.println("NUMBER");
				running = running + "" + c;
			}
			else if (isLetter(c)) {
				parsed.add(running);
				running = "";
				char next = e.charAt(i + 1);
				char two = e.charAt(i + 2);
				char three = e.charAt(i + 3);
				char four = e.charAt(i + 4);
				if (isFunction("" + c + "" + next + "" + two + "" + three + "" + four)) {
					i += 4;
					running = "" + c + "" + next + "" + two + "" + three + "" + four; 
					parsed.add(running);
					running = "";
				}
				if (isFunction("" + c + "" + next + "" + two + "" + three)) {
					i += 3;
					running = "" + c + "" + next + "" + two + "" + three;
					parsed.add(running);
					running = "";
				}
				else if (isFunction("" + c + "" + next + "" + two)) {
					i += 2;
					running = "" + c + "" + next + "" + two;
					parsed.add(running);
					running = "";
				}
				else if (isFunction("" + c + "" + next)) {
					i++;
					running = "" + c + "" + next;
					parsed.add(running);
					running = "";
				}
				
				else if (isConstant("" + c)) {
					running = "" + c;
					parsed.add(running);
					running = "";
				}
				else if (isConstant("" + c + "" + next)) {
					i++;
					running = "" + c + "" + next;
					parsed.add(running);
					running = "";
				}
				else if (isConstant("" + c + "" + next + "" + two)) {
					i += 2;
					running = "" + c + "" + next + "" + two;
					parsed.add(running);
					running = "";
				}
			}
		}
		
		parsed.add(running);
		
		for (int i = 0; i < parsed.size(); ) {
			if (parsed.get(i).equals(" ")) parsed.remove(i);
			else if (parsed.get(i).equals("")) parsed.remove(i);
			else i++;
		}
		
	//	System.out.println(parsed);
		return parsed;
	}
	
	static void help() {
		System.out.println("Type in an expression to be solved.");
	}
	
	static void asset() {
		System.out.println("SUPPORTED FUNCTIONS");
		System.out.println("sin cos tan sec csc cot exp ln log phi sigma");
	}
}


public class calc_all { // hmmm...
	static double ans;
	
	public static void main(String[] args) {
		info.init();
		
		String F = "factorize(";
		
		Scanner x = new Scanner(System.in);
		while (true) {
			String s = x.nextLine();
			s = s.toLowerCase();
		//	System.out.print(s); 
			
			if (s.equals("help")) {
				info.help();
				continue;
			}
			else if (s.length() >= 4 && s.substring(0, 4).equals("help")) {
				info.help(s.substring(5));
				continue;
			}
			else if (s.equals("functions")) {
				info.asset();
				continue;
			}
			
			// functions that return structures
			else if (s.length() >= F.length() && s.substring(0, F.length()).equals(F)) {
				if (s.charAt(s.length() - 1) != ')') s = s + ")";
				try {
					int n = Integer.parseInt(s.substring(F.length(), s.length() - 1));
					System.out.println(functions.printExpansion(n));
				}
				catch (Exception e) {
					info.error();
				}
				continue;
			}
			
			// easter eggs (nothing serious just little jokes)
			else if (info.errors.contains(s)) {
				System.out.println("You dare use my own spells against me?");
				continue;
			}
			else if (info.weird.contains(s)) {
				System.out.println("...");
				continue;
			}
			else if (s.equals("type help or functions for more information")) {
				info.idk();
				
				continue;
			}
			
			ArrayList<String> ex = parser.parse(s);
			long start = System.nanoTime();
			try {
				double eval = evaluate.evaluate(ex);
				boolean qm = ans == 4 && eval == 3;
				ans = eval;
				System.out.println(" = " + ans);
				long end = System.nanoTime();
				if (qm) System.out.println("quick maths"); // another easter egg
				System.out.println(">> " + ((end - start) / 1000000) + "ms");
			}
			catch (Exception ee) {
				info.error();
			}
				
		}
	}
}

