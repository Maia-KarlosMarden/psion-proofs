const fs=require('fs'); eval(fs.readFileSync(__dirname+'/coreT.js','utf8')+';globalThis.simulate=simulate;');
const ref=JSON.parse(fs.readFileSync(__dirname+'/../../P14TPROG_results.json','utf8'));
const H={H1:[['cool',1e-5]],H2:[['cool',0.01],['reheat',2],['cool',1e-5]],H3:[['cool',1e6],['reheat',1e11],['cool',1e-5]],
 H4:[['cool',0.05],['reheat',5],['cool',0.05],['reheat',5],['cool',1e-5]],H5:[['cool',1e-4],['reheat',0.05],['cool',1e-5]]};
let worst=0,wk='';
for(const h in H) for(const [lab,eta] of [['eta0',{}],['eta1e-9',{N8:1e-9,e10:1e-9}]]){
  const r=simulate(1e12,H[h],eta).final, p=ref[h+'_'+lab];
  for(const s of ['topo','N8','e10']) for(const q of ['X','Xbar','A']){const a=r[s][q],b=p[s][q];
    const d=(b===0&&a===0)?0:Math.abs(a/b-1); if(d>worst){worst=d;wk=`${h} ${lab} ${s} ${q} js=${a} py=${b}`;}}}
console.log('C4 maior desvio relativo JS vs Python:',worst.toExponential(2),worst<1e-6?'OK':'FALHA','|',wk);
