dtmc

const double p;
const double q;

module test

	// local state
	s : [0..8] init 0;
	
	[] s=0 -> p : (s'=7) + (1-p) : (s'=8);
    [] s=7 -> p : (s'=1) + (1-p) : (s'=2);
    [] s=8 -> p : (s'=3) + (1-p) : (s'=4);
	[] s=1 -> 0.25 : (s'=4) + 0.75 : (s'=5);
	[] s=2 -> 0.5 : (s'=4) + 0.5 : (s'=5);
	[] s=3 -> 0.2 : (s'=4) + 0.8 : (s'=5);
	[] s=4 -> q : (s'=6) + (1 - q) : (s'=5);
	[] s=5 -> 1 : (s'=5);
	[] s=6 -> 1 : (s'=6);
endmodule
label "target" = s=6;
