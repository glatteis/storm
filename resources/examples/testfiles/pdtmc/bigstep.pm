dtmc

const double p;

module test

	// local state
	s : [0..5] init 0;
	
	[] s=0 -> p*p : (s'=1) + p*(1-p) : (s'=2) + (1-p) : (s'=3);
	[] s=1 -> 0.25 : (s'=4) + 0.75 : (s'=5);
	[] s=2 -> 0.5 : (s'=4) + 0.5 : (s'=5);
	[] s=3 -> 0.2 : (s'=4) + 0.8 : (s'=5);
	[] s=4 -> 1 : (s'=4);
	[] s=5 -> 1 : (s'=5);
endmodule
label "target" = s=4;