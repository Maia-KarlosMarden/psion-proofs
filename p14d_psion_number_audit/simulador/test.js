const fs=require('fs'); eval(fs.readFileSync(__dirname+'/core.js','utf8'));
const NET=JSON.parse(fs.readFileSync(__dirname+'/net.json')), REF=JSON.parse(fs.readFileSync(__dirname+'/ref.json'));
const S=makeSim(NET); let worst=0; const t0=Date.now();
const c1=S.run({cG:1e-3,g:1,law:"L1",sym:true,u8:.5,u10:.5});
console.log("C1 simetrizado: dp",c1.dp.toExponential(2),"de",c1.de.toExponential(2));
for(const k of Object.keys(REF)){
  const [cG,g,law]=k.split('_'); const r=S.run({cG:+cG,g:+g,law,sym:false,u8:.5,u10:.5});
  const ref=REF[k].dp, ratio=r.dp/ref;
  console.log(k.padEnd(14),"JS",r.dp.toExponential(3),"Py",ref.toExponential(3),"razao",ratio.toFixed(3),"Npsi",(r.Npsi-1).toExponential(1));
}
console.log("tempo medio por corrida (ms):",((Date.now()-t0)/13).toFixed(0));
const h=S.run({cG:1e-3,g:1,law:"L2",sym:false,u8:.5,u10:.5,beta:0.2});
console.log("vies hipotetico beta=0.2: dp",h.dp.toExponential(2),"de",h.de.toExponential(2));
