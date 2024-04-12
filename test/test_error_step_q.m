% Wed 21 Feb 16:37:15 CET 2024
erel0 = 1;
 t=-[0,1,2,3];
 z=rand(4,1);
 q=0.5;
 [e,dt0] = error_step_q(t',z',q,erel0,true)
 [e,dt0] = error_step_trapezoidal(t',z',erel0,true)

