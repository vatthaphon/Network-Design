function val=genConnMatrix(N,th)
val=rand(N);
less_Th=val<th;
gr_eq_Th=~less_Th;

val(less_Th)=0;
val(gr_eq_Th)=1;
end