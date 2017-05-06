module OptimalAllocation

Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("DataFrames")
Pkg.add("PyPlot")
Pkg.add("Distributions")
	using JuMP, DataFrames, PyPlot, Distributions, Ipopt

	export data, table_NLopt, table_JuMP

function dt(k)
N1 = 50
N2 = 2
d = unifrom(0.1,1)
for i in 1:N1
  teta1[i]=rand(d)
  teta2[i,1] = teta1[i]
  teta2[i,2] = 0.5*teta1[i]
end
dict=Dict("N1"=>N1,"N2"=>N2,"Skills1"=>teta1,"Siklls2"=>teta2)
return dict
end

	function max_JuMP(k)
	m=Model(solver = IpoptSolver())
	@variable(m,c1>=0)
	@variable(m,c2>=0)
  @variable(m,y1>=0)
  @variable(m,y2>=0)
  N1=50
  N2=2
	@NLobjective(m,Max,sum((log(c1[i])-(y1[i]/dt["Skills1"][i])^2+sum(log(c2[i,j])-(y2[i,j]/dt["Skills2"][i,j])^2)*0.5 for j in 1:N2)*0.02 for i in 1:N1))
	@NLconstraint(m, sum(c1[i]-y1[i]+sum((c2[i,j]-y2[i,j])*0.5 for j in 1:N2)*0.02 for i in 1:N1)-k==0)
  @NLconstraint1(m, sum((log(c1[i])-(y1[i]/dt["Skills1"][i])^2+sum(log(c2[i,j])-(y2[i,j]/dt["Skills2"][i,j])^2)*0.5 for j in 1:N2)*0.02 for i in 1:N1)
  >= sum((log(c1[s])-(y1[s]/dt["Skills1"][i])^2+sum(log(c2[s,r])-(y2[s,r]/dt["Skills2"][i,j])^2)*0.5 for j in 1:N2 & r in 1:j)*0.02 for i in 1:N1 & s in 1:i))
	print(m)
	status=solve(m)
	v=[getobjectivevalue(m);getvalue(c1);getvalue(c2);getvalue(y1);getvalue(y2)]
  return v
	end

end
