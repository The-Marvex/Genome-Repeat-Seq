import java.io.*;
import java.util.*;

public class DNARepeats{
	public static HashSet<Integer> rabinKarpIndex = new HashSet<Integer>();
	public static HashSet<Integer> zAlgoIndex = new HashSet<Integer>();
	public static HashSet<Integer> boyesMooreIndex= new HashSet<Integer>();

	public final static int d = 256;
	static int NO_OF_CHARS = 256;
	static void rabinKarpSearch(String pat, String txt, int q) 
    { 
        int M = pat.length(); 
        int N = txt.length(); 
        int i, j; 
        int p = 0; // hash value for pattern 
        int t = 0; // hash value for txt 
        int h = 1; 
      
        // The value of h would be "pow(d, M-1)%q" 
        for (i = 0; i < M-1; i++) 
            h = (h*d)%q; 
      
        // Calculate the hash value of pattern and first 
        // window of text 
        for (i = 0; i < M; i++) 
        { 
            p = (d*p + pat.charAt(i))%q; 
            t = (d*t + txt.charAt(i))%q; 
        } 
      
        // Slide the pattern over text one by one 
        for (i = 0; i <= N - M; i++) 
        { 
      
            // Check the hash values of current window of text 
            // and pattern. If the hash values match then only 
            // check for characters on by one 
            if ( p == t ) 
            { 
                /* Check for characters one by one */
                for (j = 0; j < M; j++) 
                { 
                    if (txt.charAt(i+j) != pat.charAt(j)) 
                        break; 
                } 
      
                // if p == t and pat[0...M-1] = txt[i, i+1, ...i+M-1] 
                if (j == M){

                	rabinKarpIndex.add(i); 
                } 
                    
            } 
      
            // Calculate hash value for next window of text: Remove 
            // leading digit, add trailing digit 
            if ( i < N-M ) 
            { 
                t = (d*(t - txt.charAt(i)*h) + txt.charAt(i+M))%q; 
      
                // We might get negative value of t, converting it 
                // to positive 
                if (t < 0) 
                t = (t + q); 
            } 
        } 
    }
    static int max (int a, int b) { return (a > b)? a: b; } 
  
     //The preprocessing function for Boyer Moore's 
     //bad character heuristic 
    public static void findPattern(String t, String p)
    {
        char[] text = t.toCharArray();
        char[] pattern = p.toCharArray();
        int pos = indexOf(text, pattern);
        if (pos != -1)
            boyesMooreIndex.add(pos);
    }
    /** Function to calculate index of pattern substring **/
    public static int indexOf(char[] text, char[] pattern) 
    {
        if (pattern.length == 0) 
            return 0;
        int charTable[] = makeCharTable(pattern);
        int offsetTable[] = makeOffsetTable(pattern);
        for (int i = pattern.length - 1, j; i < text.length;) 
        {
            for (j = pattern.length - 1; pattern[j] == text[i]; --i, --j) 
                     if (j == 0) 
                    return i;
 
              // i += pattern.length - j; // For naive method
              i += Math.max(offsetTable[pattern.length - 1 - j], charTable[text[i]]);
        }
        return -1;
      }
      /** Makes the jump table based on the mismatched character information **/
      private static int[] makeCharTable(char[] pattern) 
      {
        final int ALPHABET_SIZE = 256;
        int[] table = new int[ALPHABET_SIZE];
        for (int i = 0; i < table.length; ++i) 
               table[i] = pattern.length;
        for (int i = 0; i < pattern.length - 1; ++i) 
               table[pattern[i]] = pattern.length - 1 - i;
        return table;
      }
      /** Makes the jump table based on the scan offset which mismatch occurs. **/
      private static int[] makeOffsetTable(char[] pattern) 
      {
        int[] table = new int[pattern.length];
        int lastPrefixPosition = pattern.length;
        for (int i = pattern.length - 1; i >= 0; --i) 
        {
            if (isPrefix(pattern, i + 1)) 
                   lastPrefixPosition = i + 1;
              table[pattern.length - 1 - i] = lastPrefixPosition - i + pattern.length - 1;
        }
        for (int i = 0; i < pattern.length - 1; ++i) 
        {
              int slen = suffixLength(pattern, i);
              table[slen] = pattern.length - 1 - i + slen;
        }
        return table;
    }
    /** function to check if needle[p:end] a prefix of pattern **/
    private static boolean isPrefix(char[] pattern, int p) 
    {
        for (int i = p, j = 0; i < pattern.length; ++i, ++j) 
            if (pattern[i] != pattern[j]) 
                  return false;
        return true;
    }
    /** function to returns the maximum length of the substring ends at p and is a suffix **/
    private static int suffixLength(char[] pattern, int p) 
    {
        int len = 0;
        for (int i = p, j = pattern.length - 1; i >= 0 && pattern[i] == pattern[j]; --i, --j) 
               len += 1;
        return len;
    }
    private String P;
    private String T;
    private int[] Z;
    private ArrayList<Integer> match_pos = new ArrayList<Integer>();

    public void calculateZ(String P) {
	int n = P.length();
	Z = new int[n];

	int L, R, k;
	L = R = 0;
	for (int i=1; i<n; i++) {
	    if (i >= R) {
		L = R = i;
		while(R<n && P.charAt(R-L) == P.charAt(R))
		    R++;
		Z[i] = R-L;
	    } else {
		k = i - L;
		if (Z[k] < R-i+1) {
		    Z[i] = Z[k];
		} else {
		    L = i;
		    while(R<n && P.charAt(R-L) == P.charAt(R))
			R++;
		    Z[i] = R-L;
		}
	    }
	}
	
    }

    public int[] getZ() {
		return Z;
    }
    
    public void search(String T, String P) {
	this.P = P;
	this.T = T;
	int n = P.length();

	String S = P + "$" + T;
	calculateZ(S);
	for (int i = 0; i<S.length(); i++) {
	    if (Z[i] == n) 
		match_pos.add(i-n-1);
		}
    }
    
    public void printPos() {
	for (int i=0; i<match_pos.size(); i++) {
	    System.out.println(match_pos.get(i));
	}
    }
	public static void main(String[] args)throws IOException {
        String seq;
        BufferedReader read = new BufferedReader(new InputStreamReader(System.in));
        String path = "/home/luke/Data/input7.fa";
        File file = new File(path);
        BufferedReader br = new BufferedReader(new FileReader(file));
        StringBuilder completeSeq = new StringBuilder();
        while ((seq = br.readLine()) != null){
            completeSeq.append(seq);
            completeSeq.trimToSize();
        }
        System.out.println("Sequence loaded");
        DNARepeats obj = new DNARepeats();
        seq = completeSeq.toString();
        System.out.println(seq);
        String string_pattern[] = new String[]{"CAG","CGG","CTG","GAA","GCC","GCG","CCTG","ATTCT","TGGAA","GGGGCC","CCCCGCCCCGCG"};
        for (int i = 0 ; i < string_pattern.length ; i++) {
			rabinKarpSearch(string_pattern[i], seq, 101);
			findPattern(seq, string_pattern[i]);
			obj.search(string_pattern[i], seq);
			int lastIndex = 0;
			int count = 0;
			int max_val = -1;
			List<Integer> listIndex = new ArrayList<Integer>(rabinKarpIndex);
			Collections.sort(listIndex);
			int x = 0;
			for (int val : listIndex) {
				//System.out.println(val+" "+lastIndex+" "+string_pattern[i].length());
				if(x == 0){
					lastIndex = val;
					x++;
				}
				else{
					if ((lastIndex+string_pattern[i].length())==val) {
						count++;
						if (max_val<count) {
						max_val = count;
						}
					}
					else{
						count = 0;
					}
				}
				lastIndex = val;
			}
			System.out.println("The Number of continuous repeats of "+string_pattern[i]+" are "+(max_val+1));
			//System.out.println("The Number of repeats of "+string_pattern[i]+" are "+obj.boyesMooreIndex.size());
			//System.out.println("The Number of repeats of "+string_pattern[i]+" are "+obj.zAlgoIndex.size());
			obj.rabinKarpIndex.clear();
			obj.boyesMooreIndex.clear();
			obj.zAlgoIndex.clear();			        	
        }
    }
}
