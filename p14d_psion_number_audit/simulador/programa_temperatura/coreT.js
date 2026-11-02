/*CORE-START*/
const MPL=1.2209e19, INV_B=Math.log(610/0.511)/2, LAM=0.511e-3*Math.exp(10*INV_B), ALPHA=1/137.036;
const mLev=l=>LAM*Math.exp(-l*INV_B);
function gstar(T){const c=[[80,106.75],[4,86.25],[1,75.75],[0.15,61.75],[0.1,17.25],[5e-4,10.75]];
  for(const [cut,g] of c) if(T>cut) return g; return 3.91;}
function i0f(x){const t=(x/3.75)**2; return 1+t*(3.5156229+t*(3.0899424+t*(1.2067492+t*(0.2659732+t*(0.0360768+t*0.0045813)))));}
function i1f(x){const t=(x/3.75)**2; return x*(0.5+t*(0.87890594+t*(0.51498869+t*(0.15084934+t*(0.02658733+t*(0.00301532+t*0.00032411))))));}
function k0e(x){ if(x<=2){const t=x*x/4;
  return Math.exp(x)*(-Math.log(x/2)*i0f(x)-0.57721566+t*(0.42278420+t*(0.23069756+t*(0.03488590+t*(0.00262698+t*(0.00010750+t*0.0000074))))));}
  const t=2/x; return (1.25331414+t*(-0.07832358+t*(0.02189568+t*(-0.01062446+t*(0.00587872+t*(-0.00251540+t*0.00053208))))))/Math.sqrt(x);}
function k1e(x){ if(x<=2){const t=x*x/4;
  return Math.exp(x)*(x*Math.log(x/2)*i1f(x)+1+t*(0.15443144+t*(-0.67278579+t*(-0.18156897+t*(-0.01919402+t*(-0.00110404-t*0.00004686))))))/x;}
  const t=2/x; return (1.25331414+t*(0.23498619+t*(-0.03655620+t*(0.01504268+t*(-0.00780353+t*(0.00325614-t*0.00068245))))))/Math.sqrt(x);}
function feq(x){ if(x<1e-6) return 1; if(x>700) return 0; return 0.5*x*x*(k0e(x)+2*k1e(x)/x)*Math.exp(-x); }
function makeSectors(o){ o=o||{}; const mt=mLev(o.lTop==null?0.8:o.lTop), m8=o.m8||0.938;
  return [{name:'topo',m:mt,g:2,sv:(o.kTop==null?1:o.kTop)/(mt*mt)},
          {name:'N8',m:m8,g:2,sv:1/0.135**2},
          {name:'e10',m:0.511e-3,g:2,sv:Math.PI*ALPHA*ALPHA/(0.511e-3)**2}]; }
function Yeq(sec,T){const s=2*Math.PI**2/45*gstar(T)*T**3; return sec.g*T**3/Math.PI**2*feq(sec.m/T)/s;}
// prog: [['cool',Tend],['reheat',T2],...]; eta: {N8:..., e10:...} (n/n_gamma); retorna final e trajetoria
function simulate(T0,prog,eta,opt){ opt=opt||{}; const dl=opt.dl||0.004, secs=makeSectors(opt), every=opt.every||25;
  let T=T0, lna=0, dil=1, k=0; const st=secs.map(s=>{const A=((eta&&eta[s.name])||0)/7.04, e=Yeq(s,T); return [A,(-A+Math.sqrt(A*A+4*e*e))/2];});
  const traj=[]; const rec=()=>traj.push({lna,T,v:st.map(o=>[7.04*(o[1]+o[0]),7.04*o[1]]),eq:secs.map(s=>7.04*Yeq(s,T))});
  rec();
  for(const [kind,Tt] of prog){
    if(kind==='reheat'){ if(Tt>T){const f=gstar(T)*T**3/(gstar(Tt)*Tt**3); dil*=f; st.forEach(o=>{o[0]*=f;o[1]*=f;}); T=Tt; rec();} continue; }
    while(T>Tt){ const Tn=Math.max(Tt,T*Math.exp(-dl)); lna+=Math.log(T/Tn); T=Tn;
      const H=1.66*Math.sqrt(gstar(T))*T*T/MPL, s=2*Math.PI**2/45*gstar(T)*T**3, h=dl;
      secs.forEach((sec,i)=>{const o=st[i], kk=s*sec.sv/H*h, e=Yeq(sec,T), b=1+kk*o[0], c=o[1]+kk*e*e;
        o[1]=2*c/(b+Math.sqrt(b*b+4*kk*c));});
      if((++k)%every===0) rec(); }
    rec(); }
  const fin={}; secs.forEach((sec,i)=>{const [A,yb]=st[i]; fin[sec.name]={X:7.04*(yb+A),Xbar:7.04*yb,A:7.04*A,Oh2:2.74e8*sec.m*(2*yb+A),m:sec.m};});
  return {final:fin,dil,traj,secs:secs.map(s=>({name:s.name,m:s.m}))};
}
/*CORE-END*/
