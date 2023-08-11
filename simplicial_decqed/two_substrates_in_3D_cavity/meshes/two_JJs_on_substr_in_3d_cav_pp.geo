// The parameters for the 3d cavity
L = 1000;
H = 550;
R = 80;
// Parameters for the two identical substrates
s = 50;
w = 25;
t = 3;
l = 600;
// Parameters for the JJ pads
c1 = 6;
c2 = 3;
d  = 0.5;


// The cavity
Point(1) = {-L/2, -R, -H/2};
Point(2) = {-L/2-R, 0, -H/2};
Point(3) = {-L/2, R, -H/2};
Point(4) = {L/2, R, -H/2};
Point(5) = {L/2+R, 0, -H/2};
Point(6) = {L/2, -R, -H/2};
Point(7) = {-L/2, -R, H/2};
Point(8) = {-L/2-R, 0, H/2};
Point(9) = {-L/2, R, H/2};
Point(10) = {L/2, R, H/2};
Point(11) = {L/2+R, 0, H/2};
Point(12) = {L/2, -R, H/2};
Point(13) = {-L/2, 0, -H/2};
Point(14) = {L/2, 0, -H/2};
Point(15) = {-L/2, 0, H/2};
Point(16) = {L/2, 0, H/2};

Circle(1) = {1, 13, 2};
Circle(2) = {2, 13, 3};
Line(3)   = {3, 4};
Circle(4) = {4, 14, 5};
Circle(5) = {5, 14, 6};
Line(6)   = {6, 1};
Circle(7) = {7, 15, 8};
Circle(8) = {8, 15, 9};
Line(9)   = {9, 10};
Circle(10) = {10, 16, 11};
Circle(11) = {11, 16, 12};
Line(12)   = {12, 7};
Line(13)   = {1, 7};
Line(14)   = {2, 8};
Line(15)   = {3, 9};
Line(16)   = {4, 10};
Line(17)   = {5, 11};
Line(18)   = {6, 12};

Line Loop(1) = {1,2,3,4,5,6};
Line Loop(2) = {7,8,9,10,11,12};
Line Loop(3) = {6,13,-12,-18};
Line Loop(4) = {-3,15,9,-16};
Line Loop(5) = {13,7,-14,-1};
Line Loop(6) = {14,8,-15,-2};
Line Loop(7) = {16,10,-17,-4};
Line Loop(8) = {17,11,-18,-5};


// Substrate 1
Point(17) = {-l/2-w/2,-s/2, t/2, 1};
Point(18) = {-l/2+w/2,-s/2, t/2, 1};
Point(19) = {-l/2+w/2,-s/2,-t/2, 1};
Point(20) = {-l/2-w/2,-s/2,-t/2, 1};
Point(21) = {-l/2-w/2, s/2, t/2, 1};
Point(22) = {-l/2+w/2, s/2, t/2, 1};
Point(23) = {-l/2+w/2, s/2,-t/2, 1};
Point(24) = {-l/2-w/2, s/2,-t/2, 1};

Line(19) = {17,18};
Line(20) = {18,19};
Line(21) = {19,20};
Line(22) = {20,17};
Line(23) = {21,22};
Line(24) = {22,23};
Line(25) = {23,24};
Line(26) = {24,21};
Line(27) = {17,21};
Line(28) = {20,24};
Line(29) = {18,22};
Line(30) = {19,23};

Line Loop(9)  = {19,20,21,22}; // Order matters
Line Loop(10) = {23,24,25,26};
Line Loop(11) = {19,29,-23,-27};
Line Loop(12) = {-21,30,25,-28};
Line Loop(13) = {22,27,-26,-28};
Line Loop(14) = {-20,29,24,-30};


// Substrate 2
Point(25) = {l/2-w/2,-s/2, t/2, 1};
Point(26) = {l/2+w/2,-s/2, t/2, 1};
Point(27) = {l/2+w/2,-s/2,-t/2, 1};
Point(28) = {l/2-w/2,-s/2,-t/2, 1};
Point(29) = {l/2-w/2, s/2, t/2, 1};
Point(30) = {l/2+w/2, s/2, t/2, 1};
Point(31) = {l/2+w/2, s/2,-t/2, 1};
Point(32) = {l/2-w/2, s/2,-t/2, 1};

Line(31) = {25,26};
Line(32) = {26,27};
Line(33) = {27,28};
Line(34) = {28,25};
Line(35) = {29,30};
Line(36) = {30,31};
Line(37) = {31,32};
Line(38) = {32,29};
Line(39) = {25,29};
Line(40) = {28,32};
Line(41) = {26,30};
Line(42) = {27,31};

Line Loop(15)  = {31,32,33,34}; // Order matters
Line Loop(16) = {35,36,37,38};
Line Loop(17) = {31,41,-35,-39};
Line Loop(18) = {-33,42,37,-40};
Line Loop(19) = {34,39,-38,-40};
Line Loop(20) = {-32,41,36,-42};

// JJ1 (with capacitor pads)
Point(33) = {-l/2-c1/2, -d/2-c2, t/2};
Point(34) = {-l/2-c1/2, -d/2, t/2};
Point(35) = {-l/2+c1/2, -d/2, t/2};
Point(36) = {-l/2+c1/2, -d/2-c2, t/2};
Point(37) = {-l/2-c1/2, d/2, t/2};
Point(38) = {-l/2-c1/2, d/2+c2, t/2};
Point(39) = {-l/2+c1/2, d/2+c2, t/2};
Point(40) = {-l/2+c1/2, d/2, t/2};
Point(41) = {-l/2, -d/2, t/2};
Point(42) = {-l/2, d/2, t/2};

Line(43) = {33,34};
Line(44) = {34,41};
Line(45) = {41,35};
Line(46) = {35,36};
Line(47) = {36,33};
Line(48) = {37,38};
Line(49) = {38,39};
Line(50) = {39,40};
Line(51) = {40,42};
Line(52) = {42,37};
Line(53) = {41,42};

Line Loop(21) = {43,44,45,46,47};
Line Loop(22) = {48,49,50,51,52};

// JJ2 (with capacitor pads)
Point(43) = {l/2-c1/2, -d/2-c2, t/2};
Point(44) = {l/2-c1/2, -d/2, t/2};
Point(45) = {l/2+c1/2, -d/2, t/2};
Point(46) = {l/2+c1/2, -d/2-c2, t/2};
Point(47) = {l/2-c1/2, d/2, t/2};
Point(48) = {l/2-c1/2, d/2+c2, t/2};
Point(49) = {l/2+c1/2, d/2+c2, t/2};
Point(50) = {l/2+c1/2, d/2, t/2};
Point(51) = {l/2, -d/2, t/2};
Point(52) = {l/2, d/2, t/2};

Line(54) = {43,44};
Line(55) = {44,51};
Line(56) = {51,45};
Line(57) = {45,46};
Line(58) = {46,43};
Line(59) = {47,48};
Line(60) = {48,49};
Line(61) = {49,50};
Line(62) = {50,52};
Line(63) = {52,47};
Line(64) = {51,52};

Line Loop(23) = {54,55,56,57,58};
Line Loop(24) = {59,60,61,62,63};

// Create nodes on the lines
// the cavity
Transfinite Line{1,2,4,5,7,8,10,11} = 4;
Transfinite Line{3,6,9,12} = 14;
Transfinite Line{13,14,15,16,17,18} = 9;
// the substrates
Transfinite Line{19,21,23,25,31,33,35,37} = 6;
Transfinite Line{27,28,29,30,39,40,41,42} = 10;
Transfinite Line{20,22,24,26,32,34,36,38} = 5;
// the capacitor pads
Transfinite Line{43,46,48,50,54,57,59,61} = 5;
Transfinite Line{47,49,58,60} = 7;
Transfinite Line{44,45,51,52,55,56,62,63} = 4;
// the JJ lines
Transfinite Line{53,64} = 4;

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7};
Physical Line(8) = {8};
Physical Line(9) = {9};
Physical Line(10) = {10};
Physical Line(11) = {11};
Physical Line(12) = {12};
Physical Line(13) = {13};
Physical Line(14) = {14};
Physical Line(15) = {15};
Physical Line(16) = {16};
Physical Line(17) = {17};
Physical Line(18) = {18};
Physical Line(19) = {19};
Physical Line(20) = {20};
Physical Line(21) = {21};
Physical Line(22) = {22};
Physical Line(23) = {23};
Physical Line(24) = {24};
Physical Line(25) = {25};
Physical Line(26) = {26};
Physical Line(27) = {27};
Physical Line(28) = {28};
Physical Line(29) = {29};
Physical Line(30) = {30};
Physical Line(31) = {31};
Physical Line(32) = {32};
Physical Line(33) = {33};
Physical Line(34) = {34};
Physical Line(35) = {35};
Physical Line(36) = {36};
Physical Line(37) = {37};
Physical Line(38) = {38};
Physical Line(39) = {39};
Physical Line(40) = {40};
Physical Line(41) = {41};
Physical Line(42) = {42};
Physical Line(43) = {43};
Physical Line(44) = {44};
Physical Line(45) = {45};
Physical Line(46) = {46};
Physical Line(47) = {47};
Physical Line(48) = {48};
Physical Line(49) = {49};
Physical Line(50) = {50};
Physical Line(51) = {51};
Physical Line(52) = {52};
Physical Line("JJ1") = {53};
Physical Line(54) = {54};
Physical Line(55) = {55};
Physical Line(56) = {56};
Physical Line(57) = {57};
Physical Line(58) = {58};
Physical Line(59) = {59};
Physical Line(60) = {60};
Physical Line(61) = {61};
Physical Line(62) = {62};
Physical Line(63) = {63};
Physical Line("JJ2") = {64};

Plane Surface(1) = {1};	// Create a surface out of the line loop
Plane Surface(2) = {2};	// Create a surface out of the line loop
Plane Surface(3) = {3};	// Create a surface out of the line loop
Plane Surface(4) = {4};	// Create a surface out of the line loop
Plane Surface(5) = {5};	// Create a surface out of the line loop
Plane Surface(6) = {6};	// Create a surface out of the line loop
Plane Surface(7) = {7};	// Create a surface out of the line loop
Plane Surface(8) = {8};	// Create a surface out of the line loop
Plane Surface(9) = {9};	// Create a surface out of the line loop
Plane Surface(10) = {10};	// Create a surface out of the line loop
Plane Surface(11) = {11,21,22};	// Create a surface out of the line loop
Plane Surface(12) = {12};	// Create a surface out of the line loop
Plane Surface(13) = {13};	// Create a surface out of the line loop
Plane Surface(14) = {14};	// Create a surface out of the line loop
Plane Surface(15) = {15};	// Create a surface out of the line loop
Plane Surface(16) = {16};	// Create a surface out of the line loop
Plane Surface(17) = {17,23,24};	// Create a surface out of the line loop
Plane Surface(18) = {18};	// Create a surface out of the line loop
Plane Surface(19) = {19};	// Create a surface out of the line loop
Plane Surface(20) = {20};	// Create a surface out of the line loop
Plane Surface(21) = {21};	// Create a surface out of the line loop
Plane Surface(22) = {22};	// Create a surface out of the line loop
Plane Surface(23) = {23};	// Create a surface out of the line loop
Plane Surface(24) = {24};	// Create a surface out of the line loop


Line{53} In Surface{11};
Line{64} In Surface{17};

Physical Surface("boundary1") = {1}; // Make plane surface 1 a physical surface
Physical Surface("boundary2") = {2}; // Make plane surface 1 a physical surface
Physical Surface("boundary3") = {3}; // Make plane surface 1 a physical surface
Physical Surface("boundary4") = {4}; // Make plane surface 1 a physical surface
Physical Surface("boundary5") = {5}; // Make plane surface 1 a physical surface
Physical Surface("boundary6") = {6}; // Make plane surface 1 a physical surface
Physical Surface("boundary7") = {7}; // Make plane surface 1 a physical surface
Physical Surface("boundary8") = {8}; // Make plane surface 1 a physical surface
Physical Surface(9) = {9}; // Make plane surface 1 a physical surface
Physical Surface(10) = {10}; // Make plane surface 1 a physical surface
Physical Surface("substrate 1's surface") = {11}; // Make plane surface 1 a physical surface
Physical Surface(12) = {12}; // Make plane surface 1 a physical surface
Physical Surface(13) = {13}; // Make plane surface 1 a physical surface
Physical Surface(14) = {14}; // Make plane surface 1 a physical surface
Physical Surface(15) = {15}; // Make plane surface 1 a physical surface
Physical Surface(16) = {16}; // Make plane surface 1 a physical surface
Physical Surface("substrate 2's surface") = {17}; // Make plane surface 1 a physical surface
Physical Surface(18) = {18}; // Make plane surface 1 a physical surface
Physical Surface(19) = {19}; // Make plane surface 1 a physical surface
Physical Surface(20) = {20}; // Make plane surface 1 a physical surface
Physical Surface("superconductor1") = {21}; // Make plane surface 1 a physical surface
Physical Surface("superconductor2") = {22}; // Make plane surface 1 a physical surface
Physical Surface("superconductor3") = {23}; // Make plane surface 1 a physical surface
Physical Surface("superconductor4") = {24}; // Make plane surface 1 a physical surface



Surface Loop(1) = {1,-3,5,6,4,7,8,-2}; 
Surface Loop(2) = {-12,13,10,-14,-9,11,-21,-22}; 
Surface Loop(3) = {-18,19,16,-20,-15,17,-23,-24};
Volume(1) = {1,2,3};
Volume(2) = {2};
Volume(3) = {3};
Physical Volume ("vacuum") = {1};
Physical Volume ("substrate1") = {2};
Physical Volume ("substrate2") = {3};


