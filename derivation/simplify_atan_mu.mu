//Mon 27 Aug 19:44:48 CEST 2018
simplify_atan_mu := proc(x)
local _rewrite;
begin
	_rewrite := proc(y)
	local r;
	begin
		r := Rule(atan(`#X`) + atan(`#Y`), atan2( (#X + #Y),(1 - #X*#Y) ) );
		tmp := Rule::apply(r,y);
		if tmp <> FAIL then
       		     y:= tmp;
         	end;
		return(y);
	end;
    	misc::maprec(x, TRUE = _rewrite);
end_proc:

