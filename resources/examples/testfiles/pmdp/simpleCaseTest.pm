mdp

// parameters
const double p;
const double q;

module test

	// local state
	s : [0..6] init 0;

  // s0 alpha
	[] s=0 -> p : (s'=1) + 1-p : (s'=2);
  // s0 beta
	[] s=0 -> q : (s'=3) + 1-q : (s'=4);

	// s1 to s4
	[]s=1 -> 0.9 : (s'=5) + 0.1 : (s'=6);
  []s=2 -> 0.8 : (s'=1) + 0.2 : (s'=6);
  []s=3 -> 0.2 : (s'=2) + 0.8 : (s'=6);
  []s=4 -> 0.1 : (s'=3) + 0.9 : (s'=6);

  // abs states s5 and s6
	[] s=5 -> 1 : (s'=5);
	[] s=6 -> 1 : (s'=6);


endmodule
