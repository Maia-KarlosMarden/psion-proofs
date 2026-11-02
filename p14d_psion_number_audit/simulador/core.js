/*CORE-START*/
function makeSim(NET){
  const NT=20, NS=6+2*NT;
  const TR=NET.TR, CJ=NET.CJ, RE=NET.REARR;
  const iP=TR.findIndex(t=>t.join()==="0,4,5"), iPb=TR.findIndex(t=>t.join()==="0,1,2");
  function weights(g,law,sym){
    const H=NET.H0.map((r,i)=>r.map((x,j)=>x+g*NET.O[i][j]));
    const d=H.map((r,i)=>r[i]);
    let kb=law==="L1"?d.map(x=>x*x):(()=>{const m=Math.min(...d);return d.map(x=>Math.exp(-(x-m)));})();
    const mean=kb.reduce((a,b)=>a+b,0)/NT; kb=kb.map(x=>x/mean);
    let W=H.map((r,i)=>r.map((x,j)=>i===j?0:x*x));
    if(sym){
      kb=kb.map((x,i)=>0.5*(x+kb[CJ[i]]));
      W=W.map((r,i)=>r.map((x,j)=>0.5*(x+W[CJ[i]][CJ[j]])));
    }
    return {kb,W};
  }
  function gammaH(T,cG){return cG*NET.G4*T**3*NET.MPL/(1.66*Math.sqrt(NET.GSTAR)*NET.LAM**4);}
  function run(p){
    const {kb,W}=weights(p.g,p.law,p.sym);
    const u=[p.u8,p.u10];
    const bias=[new Array(NT).fill(1),new Array(NT).fill(1)];
    if(p.beta){bias[0][iP]=1+p.beta;bias[0][iPb]=1-p.beta;bias[1][iPb]=1+p.beta;bias[1][iP]=1-p.beta;}
    const Wrow=W.map(r=>r.reduce((a,b)=>a+b,0));
    const K=NET.K_EQ;
    let Y=new Float64Array(NS); for(let k=0;k<6;k++)Y[k]=1/6;
    const xend=Math.log(NET.T_RH/NET.T_END), nst=p.steps||1800, dx=xend/nst;
    const F=new Float64Array(NS), Jm=[]; for(let i=0;i<NS;i++)Jm.push(new Float64Array(NS));
    const traj=[];
    function evalFJ(Yv,r){
      F.fill(0); for(let i=0;i<NS;i++)Jm[i].fill(0);
      for(let l=0;l<2;l++){
        const off=6+l*NT;
        for(let i=0;i<NT;i++){
          const [a,b,c]=TR[i], k=kb[i]*u[l]*bias[l][i], t=off+i;
          const J=k*(Yv[a]*Yv[b]*Yv[c]-Yv[t]/K);
          F[t]+=J; F[a]-=J; F[b]-=J; F[c]-=J;
          const dA=k*Yv[b]*Yv[c], dB=k*Yv[a]*Yv[c], dC=k*Yv[a]*Yv[b], dT=-k/K;
          for(const [row,s] of [[t,1],[a,-1],[b,-1],[c,-1]]){
            Jm[row][a]+=s*dA; Jm[row][b]+=s*dB; Jm[row][c]+=s*dC; Jm[row][t]+=s*dT;
          }
          for(let jj=0;jj<NT;jj++){ if(W[i][jj]===0)continue;
            F[t]+=W[i][jj]*(Yv[off+jj]-Yv[t]); Jm[t][off+jj]+=W[i][jj]; }
          Jm[t][t]-=Wrow[i];
        }
        for(const [a,b,c,d] of RE){
          const A=off+a,B=off+b,C=off+c,D=off+d;
          const Jr=Yv[A]*Yv[B]-Yv[C]*Yv[D];
          F[A]-=Jr;F[B]-=Jr;F[C]+=Jr;F[D]+=Jr;
          for(const [row,s] of [[A,-1],[B,-1],[C,1],[D,1]]){
            Jm[row][A]+=s*Yv[B]; Jm[row][B]+=s*Yv[A]; Jm[row][C]-=s*Yv[D]; Jm[row][D]-=s*Yv[C];
          }
        }
      }
      for(let i=0;i<NS;i++){F[i]*=r; for(let j=0;j<NS;j++)Jm[i][j]*=r;}
    }
    function solve(M,v){ // Gauss com pivoteamento, in place
      const n=v.length;
      for(let c=0;c<n;c++){
        let p=c,mx=Math.abs(M[c][c]); for(let r=c+1;r<n;r++){const a=Math.abs(M[r][c]); if(a>mx){mx=a;p=r;}}
        if(p!==c){[M[p],M[c]]=[M[c],M[p]];[v[p],v[c]]=[v[c],v[p]];}
        const pv=M[c][c];
        for(let r=c+1;r<n;r++){const f=M[r][c]/pv; if(f===0)continue; const Mr=M[r],Mc=M[c]; for(let k=c;k<n;k++)Mr[k]-=f*Mc[k]; v[r]-=f*v[c];}
      }
      for(let c=n-1;c>=0;c--){let s=v[c]; for(let k=c+1;k<n;k++)s-=M[c][k]*v[k]; v[c]=s/M[c][c];}
      return v;
    }
    const record=(x)=>{const off8=6,off10=6+NT;
      traj.push({T:NET.T_RH*Math.exp(-x),p:Y[off8+iP],pb:Y[off8+iPb],e:Y[off10+iPb],eb:Y[off10+iP],
        dp:Y[off8+iP]-Y[off8+iPb],de:Y[off10+iPb]-Y[off10+iP]});};
    record(0);
    for(let n=1;n<=nst;n++){
      const x=n*dx, r=gammaH(NET.T_RH*Math.exp(-x),p.cG);
      const Yn=Float64Array.from(Y);
      for(let it=0;it<4;it++){
        evalFJ(Y,r);
        const G=new Float64Array(NS), M=[];
        for(let i=0;i<NS;i++){G[i]=-(Y[i]-Yn[i]-dx*F[i]); const row=new Float64Array(NS); for(let j=0;j<NS;j++)row[j]=-dx*Jm[i][j]; row[i]+=1; M.push(row);}
        const dY=solve(M,G); let mx=0;
        for(let i=0;i<NS;i++){Y[i]+=dY[i]; mx=Math.max(mx,Math.abs(dY[i]));}
        if(mx<1e-15)break;
      }
      if(n%12===0||n===nst)record(x);
    }
    let Npsi=0; for(let k=0;k<6;k++)Npsi+=Y[k]; for(let i=6;i<NS;i++)Npsi+=3*Y[i];
    let Q=0; for(let k=0;k<6;k++)Q+=NET.Q9[k]/9*Y[k];
    for(let l=0;l<2;l++)for(let i=0;i<NT;i++){Q+=TR[i].reduce((a,k)=>a+NET.Q9[k],0)/9*Y[6+l*NT+i];}
    const last=traj[traj.length-1];
    return {traj,dp:last.dp,de:last.de,eta:7.04*NET.Y0*last.dp,Npsi,Q};
  }
  return {run};
}
/*CORE-END*/
