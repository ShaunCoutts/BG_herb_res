## Dynamic program calling script that pulls in all functions and sets up the runs 

# DEFINE ACTION STATE INDEX CONST to make the code a bit more readable 
const HERB = 1;
const SPOT = 2;
const CROP = 3;
const BANK = 4;
# sub-actions index, i.e. within each action
const NO_H = 1;
const ONE_H = 2;
const MUL_H = 3;
const NO_SP = 1;
const YS_SP = 2;
const WHEAT = 1;
const ALT = 2;
const FAL = 3;
const NO_SB = 1;
const YS_SB = 2;
# DEFINE STATE INDEX to make code more readable
const S_g = 1;
const S_gRR = 1;
const S_gRr = 2;
const S_grr = 3;
const S_N = 2;
const S_pRR = 3;
const S_pRr = 4;
# DOMINANCE OF TSR
const resist_G = ["RR", "Rr"];

# parameters 
seed_sur_base = 0.5;
seed_sur_act = 0.4;



# ACTION SPACE 
sub_actions = [[NO_H, ONE_H, MUL_H], 
  [NO_SP, YS_SP], 
  [WHEAT, ALT, FAL], 
  [NO_SB, YS_SB]];
  
action_space = build_action_space();

# action effects 
action_effect = ([0, 1, 2], [false, true], [WHEAT, ALT], [seed_sur_base, seed_sur_act])
