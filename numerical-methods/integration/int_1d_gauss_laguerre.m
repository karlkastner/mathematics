% Sat  4 Mar 13:12:59 CET 2023
function [w,x] = int_1d_gauss_laguerre(n)
if (n > 20)
	x=laguerre_roots(n);
	w = x./((n+1)*(laguerreL(n+1,x))).^2;
else
w = [   1.000000000000000   0.146446609406726   0.010389256501586   0.000539294705561   0.000023369972386   0.000000898547906   0.000000031703155   0.000000001048001   0.000000000032909   0.000000000000991   0.000000000000029   0.000000000000001   0.000000000000000   0.000000000000000   0.000000000000000   0.000000000000000   0.000000000000000   0.000000000000000   0.000000000000000   0.000000000000000
                   0   0.853553390593274   0.278517733569238   0.038887908515005   0.003611758679922   0.000261017202815   0.000015865464349   0.000000848574672   0.000000041107693   0.000000001839565   0.000000000077126   0.000000000003062   0.000000000000116   0.000000000000004   0.000000000000000   0.000000000000000   0.000000000000000   0.000000000000000   0.000000000000000   0.000000000000000
                   0                   0   0.711093009929174   0.357418692437795   0.075942449681705   0.010399197453148   0.001074010143281   0.000090765087733   0.000006592123026   0.000000424931398   0.000000024863537   0.000000001342391   0.000000000067708   0.000000000003221   0.000000000000146   0.000000000000006   0.000000000000000   0.000000000000000   0.000000000000000   0.000000000000000
                   0                   0                   0   0.603154104341639   0.398666811083170   0.113373382074045   0.020633514468720   0.002794536235228   0.000305249767093   0.000028259233496   0.000002292403880   0.000000166849388   0.000000011088416   0.000000000681931   0.000000000039219   0.000000000002127   0.000000000000110   0.000000000000005   0.000000000000000   0.000000000000000
                   0                   0                   0                   0   0.521755610582808   0.417000830772132   0.147126348657497   0.033343492261208   0.005599626610804   0.000753008388585   0.000085131224340   0.000008365055857   0.000000731731162   0.000000058015440   0.000000004227430   0.000000000286235   0.000000000018169   0.000000000001089   0.000000000000062   0.000000000000003
                   0                   0                   0                   0                   0   0.458964673949955   0.421831277861697   0.175794986637191   0.047460562765595   0.009501516975164   0.001518880846641   0.000203231592669   0.000023515473981   0.000002409585768   0.000000222631691   0.000000018810248   0.000000001469731   0.000000000107171   0.000000000007348   0.000000000000477
                   0                   0                   0                   0                   0                   0   0.409318951701294   0.418786780814322   0.199287525370956   0.062087456098770   0.014393282766806   0.002663973541822   0.000411881770449   0.000054907194631   0.000006459926754   0.000000682831935   0.000000065762729   0.000000005836095   0.000000000481668   0.000000000037255
                   0                   0                   0                   0                   0                   0                   0   0.369188589341642   0.411213980423970   0.218068287611720   0.076564453546848   0.020102381154654   0.004220396040501   0.000739890377929   0.000111674392460   0.000014844586841   0.000001768515055   0.000000191466979   0.000000019049926   0.000000001757981
                   0                   0                   0                   0                   0                   0                   0                   0   0.336126421797967   0.401119929155304   0.232781831848726   0.090449222211964   0.026432754414742   0.006192869437669   0.001212436146640   0.000204271915523   0.000030275517808   0.000004015307947   0.000000483056460   0.000000053301222
                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.308441115765015   0.389720889527879   0.244082011319567   0.103470758025313   0.033192092156752   0.008563877804723   0.001849070942902   0.000343679727190   0.000056169650353   0.000008207307158   0.000001086486081
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.284933212894201   0.377759275873181   0.252562420056996   0.115482893559796   0.040206864920209   0.011299900081184   0.002662824736384   0.000540622786214   0.000096524712061   0.000015401443820
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.264731371055439   0.365688822900673   0.258734610244263   0.126425818105831   0.047328928693631   0.014357297762471   0.003660179772038   0.000804912768018   0.000155741758652
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.247188708429963   0.353784691597779   0.263027577942101   0.136296934296371   0.054436943250556   0.017687213073684   0.004841627330753   0.001144962419886
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.231815577144853   0.342210177922749   0.265795777644318   0.145129854356211   0.061434917866635   0.021239307608808   0.006202550806425
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.218234885940092   0.331057854950872   0.267329726357356   0.152979747465585   0.068249379984170   0.024964417334436
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.206151714957795   0.320375357274591   0.267866567149199   0.159913372130890   0.074826064654208
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.195332205251761   0.310181766370141   0.267599547039336   0.166002453276131
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.185588603146912   0.300478143607192   0.266686102866043
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.176768474915886   0.291254362005948
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.168746801851140
];
x =[   1.000000000000000   3.414213562373095   6.289945082937476   9.395070912301120  12.640800844275811  15.982873980601699  19.395727862262582  22.863131736889429  26.374071890927347  29.920697012275291  33.497192847180280  37.099121044472810  40.723008669253112  44.366081711141689  48.026085572608650  51.701160339456777  55.389751789568393  59.090546435884917  62.802423154824773  66.524416522888075
                   0   0.585786437626905   2.294280360279045   4.536620296921130   7.085810005858809   9.837467418382539  12.734180291797712  15.740678641277523  18.833597788991732  21.996585811974331  25.217709339646223  28.487967250962843  31.800386302018140  35.149443660384634  38.530683306975369  41.940452648097292  45.375716536426445  48.833922716706191  52.312902441267340  55.810795779951057
                   0                   0   0.415774556783479   1.745761101158349   3.596425771040735   5.775143569104537   8.182153444562976  10.758516010181593  13.466236911092070  16.279257831388076  19.178857403296135  22.151090379424126  25.185263864530942  28.272981723921397  31.407519168576876  34.583398701469086  37.796093829188692  41.041816769536901  44.317362874014854  47.619993920858413
                   0                   0                   0   0.322547689619392   1.413403059106519   2.992736326059313   4.900353084526436   7.045905402393181   9.372985251687959  11.843785837893597  14.431613757956880  17.116855187442837  19.884635664029858  22.723381627108910  25.623894228079802  28.578729743971511  31.581771695004580  34.627927074405093  37.712905587708313  40.833057342263537
                   0                   0                   0                   0   0.263560319718141   1.188932101672621   2.567876744950756   4.266700170287718   6.204956777876078   8.330152746765746  10.605950999625996  13.006054993321017  15.510762037628396  18.104892221408303  20.776478898840001  23.515905692594721  26.315317796766692  29.168208646700784  32.069122631831164  35.013433854746694
                   0                   0                   0                   0                   0   0.222846604179261   1.026664895339196   2.251086629866116   3.783473973331470   5.552496140064274   7.509887863772992   9.621316842445554  11.861403588816085  14.210805010420266  16.654407708054475  19.180156858413341  21.778268265762676  24.440681380940497  27.160668849867875  29.932554960085223
                   0                   0                   0                   0                   0                   0   0.193043676560362   0.903701776799383   2.005135155619310   3.401433697854659   5.029284401588105   6.844525453119666   8.815001941206026  10.916499507623996  13.130282482717446  15.441527367346561  17.838284728786316  20.310767604633075  22.850850206727170  25.451702599646193
                   0                   0                   0                   0                   0                   0                   0   0.170279632305101   0.807220022742257   1.808342901740351   3.091138143034111   4.599227639418166   6.292256271126480   8.140240141541812  10.120228567680369  12.214223369682914  14.407823037146530  16.689306296498767  19.048992890727941  21.478788358695372
                   0                   0                   0                   0                   0                   0                   0                   0   0.152322227731808   0.729454549503168   1.647150545872251   2.833751337743152   4.238845929021902   5.825536218280209   7.565916226728692   9.438314336090775  11.425529319695165  13.513656201769232  15.691278516991696  17.948895410988566
                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.137793470540493   0.665418255839226   1.512610269776496   2.616597108405353   3.932102822303781   5.425336627390854   7.070338535120583   8.846685511145015  10.737990048580993  12.730881388277394  14.814293537710773
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.125796442187968   0.611757484515128   1.398564336451156   2.430801078728760   3.667622721753768   5.078018614538271   6.637829205307184   8.327825156935885  10.132423740690534  12.038802490640716
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.115722117358021   0.566131899040393   1.300629121251698   2.269949526203801   3.437086633894294   4.773513513734307   6.256725073751477   7.868618910161866   9.594392891180910
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.107142388472252   0.526857648851891   1.215595412070887   2.129283645098384   3.234256124038587   4.504205538938685   5.918141562283298   7.459017448362532
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.099747507032598   0.492691740301890   1.141057774831213   2.005193531650238   3.054353113195036   4.264215539616115   5.615174971703933
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.093307812017282   0.462696328915081   1.075176577511407   1.894888509970594   2.893651381865065   4.048925313751387
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.087649410478928   0.436150323558708   1.016520179623472   1.796175582069679   2.749199255322502
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.082638214708948   0.412490085259132   0.963957343997851   1.707306531026682
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.078169166669706   0.391268613319997   0.916582102483354
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.074158783757205   0.372126818001615
                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0   0.070539889691988

];
w = w(1:n,n);
x = x(1:n,n);
end
end