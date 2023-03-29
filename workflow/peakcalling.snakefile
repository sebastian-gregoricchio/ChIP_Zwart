#####################################
## Snakefile for ChIP peak calling ##
#####################################

from typing import List
import pathlib
import re
import numpy
import os
import pandas as pd
import math
from itertools import combinations


### Define the genomes (ZWART source)
genomes_table = pd.DataFrame({'genome_id': ['hg38', 'hg19', 'hg38_ucsc', 'hg19_ucsc', 'mm10', 'mm9', 'rn6'],
                              'fasta': ['/shared/data/Zwartlab/snakepipes_indices/hg38/BWAIndex/genome.fa',
                                        '/shared/data/Zwartlab/snakepipes_indices/hg19/BWAIndex/genome.fa',
                                        '/home/s.gregoricchio/annotations/genomes/Hg19_UCSC/one_file_fasta/hg19.fa.gz',
                                        'hg38_ucsc',
                                        '/home/s.gregoricchio/annotations/genomes/Mm10/one_file_fasta/Mus_musculus_GRCm38_Mm10_GCA_000001635.2_UCSC.fa',
                                        'mm10_ucsc',
                                        '/shared/data/Zwartlab/snakepipes_indices/Rnor_6.0/BWAIndex/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa'],
                              'blacklist': ['/shared/data/Zwartlab/snakepipes_indices/hg38/annotation/blacklist.bed',
                                            '/shared/data/Zwartlab/snakepipes_indices/hg19/annotation/blacklist.bed',
                                            '/home/s.gregoricchio/annotations/blacklist/hg38-blacklist.v2.bed',
                                            '/home/s.gregoricchio/annotations/blacklist/hg19-blacklist.v2.bed',
                                            '/home/s.gregoricchio/annotations/blacklist/mm10-blacklist.v2.bed',
                                            '/home/s.gregoricchio/annotations/blacklist/mm9-blacklist.v2.bed',
                                            '/shared/data/Zwartlab/snakepipes_indices/Rnor_6.0/annotation/Rnor_6.0.blacklist.bed'],
                              'effective_genomeSize': [2945849067, 2861343702, 2900338458, 2900338458, 2652783500, 2620345972, 2870184193],
                              'ignore_for_normalization': ["KI270728.1 KI270727.1 KI270442.1 KI270729.1 GL000225.1 KI270743.1 GL000008.2 GL000009.2 KI270747.1 KI270722.1 GL000194.1 KI270742.1 GL000205.2 GL000195.1 KI270736.1 KI270733.1 GL000224.1 GL000219.1 KI270719.1 GL000216.2 KI270712.1 KI270706.1 KI270725.1 KI270744.1 KI270734.1 GL000213.1 GL000220.1 KI270715.1 GL000218.1 KI270749.1 KI270741.1 GL000221.1 KI270716.1 KI270731.1 KI270751.1 KI270750.1 KI270519.1 GL000214.1 KI270708.1 KI270730.1 KI270438.1 KI270737.1 KI270721.1 KI270738.1 KI270748.1 KI270435.1 GL000208.1 KI270538.1 KI270756.1 KI270739.1 KI270757.1 KI270709.1 KI270746.1 KI270753.1 KI270589.1 KI270726.1 KI270735.1 KI270711.1 KI270745.1 KI270714.1 KI270732.1 KI270713.1 KI270754.1 KI270710.1 KI270717.1 KI270724.1 KI270720.1 KI270723.1 KI270718.1 KI270317.1 KI270740.1 KI270755.1 KI270707.1 KI270579.1 KI270752.1 KI270512.1 KI270322.1 GL000226.1 KI270311.1 KI270366.1 KI270511.1 KI270448.1 KI270521.1 KI270581.1 KI270582.1 KI270515.1 KI270588.1 KI270591.1 KI270522.1 KI270507.1 KI270590.1 KI270584.1 KI270320.1 KI270382.1 KI270468.1 KI270467.1 KI270362.1 KI270517.1 KI270593.1 KI270528.1 KI270587.1 KI270364.1 KI270371.1 KI270333.1 KI270374.1 KI270411.1 KI270414.1 KI270510.1 KI270390.1 KI270375.1 KI270420.1 KI270509.1 KI270315.1 KI270302.1 KI270518.1 KI270530.1 KI270304.1 KI270418.1 KI270424.1 KI270417.1 KI270508.1 KI270303.1 KI270381.1 KI270529.1 KI270425.1 KI270396.1 KI270363.1 KI270386.1 KI270465.1 KI270383.1 KI270384.1 KI270330.1 KI270372.1 KI270548.1 KI270580.1 KI270387.1 KI270391.1 KI270305.1 KI270373.1 KI270422.1 KI270316.1 KI270338.1 KI270340.1 KI270583.1 KI270334.1 KI270429.1 KI270393.1 KI270516.1 KI270389.1 KI270466.1 KI270388.1 KI270544.1 KI270310.1 KI270412.1 KI270395.1 KI270376.1 KI270337.1 KI270335.1 KI270378.1 KI270379.1 KI270329.1 KI270419.1 KI270336.1 KI270312.1 KI270539.1 KI270385.1 KI270423.1 KI270392.1 KI270394.1 X Y MT",
                                                           "X Y MT",
                                                           "KI270728.1 KI270727.1 KI270442.1 KI270729.1 GL000225.1 KI270743.1 GL000008.2 GL000009.2 KI270747.1 KI270722.1 GL000194.1 KI270742.1 GL000205.2 GL000195.1 KI270736.1 KI270733.1 GL000224.1 GL000219.1 KI270719.1 GL000216.2 KI270712.1 KI270706.1 KI270725.1 KI270744.1 KI270734.1 GL000213.1 GL000220.1 KI270715.1 GL000218.1 KI270749.1 KI270741.1 GL000221.1 KI270716.1 KI270731.1 KI270751.1 KI270750.1 KI270519.1 GL000214.1 KI270708.1 KI270730.1 KI270438.1 KI270737.1 KI270721.1 KI270738.1 KI270748.1 KI270435.1 GL000208.1 KI270538.1 KI270756.1 KI270739.1 KI270757.1 KI270709.1 KI270746.1 KI270753.1 KI270589.1 KI270726.1 KI270735.1 KI270711.1 KI270745.1 KI270714.1 KI270732.1 KI270713.1 KI270754.1 KI270710.1 KI270717.1 KI270724.1 KI270720.1 KI270723.1 KI270718.1 KI270317.1 KI270740.1 KI270755.1 KI270707.1 KI270579.1 KI270752.1 KI270512.1 KI270322.1 GL000226.1 KI270311.1 KI270366.1 KI270511.1 KI270448.1 KI270521.1 KI270581.1 KI270582.1 KI270515.1 KI270588.1 KI270591.1 KI270522.1 KI270507.1 KI270590.1 KI270584.1 KI270320.1 KI270382.1 KI270468.1 KI270467.1 KI270362.1 KI270517.1 KI270593.1 KI270528.1 KI270587.1 KI270364.1 KI270371.1 KI270333.1 KI270374.1 KI270411.1 KI270414.1 KI270510.1 KI270390.1 KI270375.1 KI270420.1 KI270509.1 KI270315.1 KI270302.1 KI270518.1 KI270530.1 KI270304.1 KI270418.1 KI270424.1 KI270417.1 KI270508.1 KI270303.1 KI270381.1 KI270529.1 KI270425.1 KI270396.1 KI270363.1 KI270386.1 KI270465.1 KI270383.1 KI270384.1 KI270330.1 KI270372.1 KI270548.1 KI270580.1 KI270387.1 KI270391.1 KI270305.1 KI270373.1 KI270422.1 KI270316.1 KI270338.1 KI270340.1 KI270583.1 KI270334.1 KI270429.1 KI270393.1 KI270516.1 KI270389.1 KI270466.1 KI270388.1 KI270544.1 KI270310.1 KI270412.1 KI270395.1 KI270376.1 KI270337.1 KI270335.1 KI270378.1 KI270379.1 KI270329.1 KI270419.1 KI270336.1 KI270312.1 KI270539.1 KI270385.1 KI270423.1 KI270392.1 KI270394.1 X Y MT M",
                                                           "X Y MT M GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605 hs37d5",
                                                           "MT M X Y JH584299.1 GL456233.1 JH584301.1 GL456211.1 GL456350.1 JH584293.1 GL456221.1 JH584297.1 JH584296.1 GL456354.1 JH584294.1 JH584298.1 JH584300.1 GL456219.1 GL456210.1 JH584303.1 JH584302.1 GL456212.1 JH584304.1 GL456379.1 GL456216.1 GL456393.1 GL456366.1 GL456367.1 GL456239.1 GL456213.1 GL456383.1 GL456385.1 GL456360.1 GL456378.1 GL456389.1 GL456372.1 GL456370.1 GL456381.1 GL456387.1 GL456390.1 GL456394.1 GL456392.1 GL456382.1 GL456359.1 GL456396.1 GL456368.1 JH584292.1 JH584295.1",
                                                           "Y X MT M NT_166325 NT_166464 NT_166452 NT_166480 NT_166448 NT_166458 NT_166443 NT_166466 NT_166476 NT_166479 NT_166478 NT_166474 NT_166471 NT_166445 NT_166465 NT_166457 NT_166470 NT_166454 NT_166472 NT_166449 NT_166481 NT_166337 NT_166459 NT_166456 NT_166473 NT_166461 NT_166475 NT_166462 NT_166444 NT_166453 NT_166446 NT_166469 NT_072868 NT_166335 NT_166467 NT_166283 NT_166338 NT_166340 NT_166442 NT_166334 NT_166286 NT_166451 NT_166336 NT_166339 NT_166290 NT_053651 NT_166450 NT_166447 NT_166468 NT_166460 NT_166477 NT_166455 NT_166291 NT_166463 NT_166433 NT_166402 NT_166327 NT_166308 NT_166309 NT_109319 NT_166282 NT_166314 NT_166303 NT_112000 NT_110857 NT_166280 NT_166375 NT_166311 NT_166307 NT_166310 NT_166323 NT_166437 NT_166374 NT_166364 NT_166439 NT_166328 NT_166438 NT_166389 NT_162750 NT_166436 NT_166372 NT_166440 NT_166326 NT_166342 NT_166333 NT_166435 NT_166434 NT_166341 NT_166376 NT_166387 NT_166281 NT_166313 NT_166380 NT_166360 NT_166441 NT_166359 NT_166386 NT_166356 NT_166357 NT_166423 NT_166384 NT_161879 NT_161928 NT_166388 NT_161919 NT_166381 NT_166367 NT_166392 NT_166406 NT_166365 NT_166379 NT_166358 NT_161913 NT_166378 NT_166382 NT_161926 NT_166345 NT_166385 NT_165789 NT_166368 NT_166405 NT_166390 NT_166373 NT_166361 NT_166348 NT_166369 NT_161898 NT_166417 NT_166410 NT_166383 NT_166362 NT_165754 NT_166366 NT_166363 NT_161868 NT_166407 NT_165793 NT_166352 NT_161925 NT_166412 NT_165792 NT_161924 NT_166422 NT_165795 NT_166354 NT_166350 NT_165796 NT_161904 NT_166370 NT_165798 NT_165791 NT_161885 NT_166424 NT_166346 NT_165794 NT_166377 NT_166418 NT_161877 NT_166351 NT_166408 NT_166349 NT_161906 NT_166391 NT_161892 NT_166415 NT_165790 NT_166420 NT_166353 NT_166344 NT_166371 NT_161895 NT_166404 NT_166413 NT_166419 NT_161916 NT_166347 NT_161875 NT_161911 NT_161897 NT_161866 NT_166409 NT_161872 NT_166403 NT_161902 NT_166414 NT_166416 NT_166421 NT_161923 NT_161937",
                                                           "X Y MT KL568162.1 KL568139.1 KL568161.1 KL568148.1 KL568157.1 KL568160.1 KL568151.1 KL568149.1 KL568141.1 KL568156.1 KL568159.1 KL568155.1 KL568140.1 KL568158.1 KL568152.1 KL568164.1 KL568150.1 KL568163.1 KL568473.1 KL568147.1 KL568409.1 KL568165.1 KL568143.1 KL568374.1 KL568142.1 KL568418.1 KL568451.1 KL568353.1 KL568153.1 KL567939.1 KL568154.1 KL568144.1 KL568410.1 KL568166.1 KL568145.1 KL568458.1 KL568432.1 KL568465.1 AABR07024382.1 AABR07024188.1 KL568405.1 KL568487.1 KL568449.1 KL567959.1 KL568438.1 AABR07024291.1 KL568414.1 AABR07024047.1 KL568318.1 KL568447.1 KL568446.1 KL568146.1 KL568319.1 KL568468.1 KL568413.1 KL568354.1 KL568362.1 KL568415.1 KL568421.1 KL568003.1 KL568434.1 KL567919.1 AABR07024421.1 KL568482.1 KL568460.1 AABR07024203.1 KL568001.1 KL568315.1 KL568336.1 KL568366.1 KL568334.1 KL568379.1 KL568346.1 AABR07024124.1 KL568419.1 KL568450.1 KL568328.1 KL568477.1 KL568470.1 KL568360.1 KL568377.1 KL568463.1 KL567940.1 KL568067.1 KL568454.1 KL568397.1 KL568125.1 KL568388.1 KL568337.1 KL568412.1 KL568126.1 KL567994.1 KL568398.1 AABR07024145.1 AABR07024104.1 KL567946.1 KL568324.1 KL568464.1 AABR07024039.1 KL568381.1 KL568440.1 AABR07024031.1 KL568417.1 KL567908.1 KL568443.1 KL568051.1 KL568101.1 KL567926.1 KL568423.1 KL568097.1 KL567914.1 KL568453.1 KL567947.1 KL568395.1 KL568445.1 KL568430.1 KL568325.1 AABR07024102.1 KL567941.1 KL567907.1 KL567891.1 KL568102.1 AABR07024120.1 KL568424.1 KL568355.1 KL568433.1 KL568080.1 KL568114.1 KL568282.1 KL568005.1 KL568019.1 KL568137.1 KL568462.1 KL568344.1 KL567901.1 KL568467.1 KL568408.1 KL568342.1 KL568079.1 KL568370.1 KL568096.1 KL568472.1 KL568431.1 KL568108.1 KL568350.1 KL568466.1 AABR07024398.1 AABR07024066.1 KL568461.1 KL568118.1 KL568357.1 AABR07024106.1 KL568480.1 KL568136.1 KL568045.1 KL567951.1 KL568116.1 KL568364.1 KL567984.1 KL568064.1 KL568292.1 KL568479.1 KL567981.1 KL568361.1 KL568471.1 KL568510.1 KL567997.1 KL568402.1 KL568452.1 KL567988.1 KL568483.1 KL568345.1 AABR07024042.1 KL568442.1 KL568057.1 KL568029.1 KL568257.1 KL568348.1 KL568187.1 KL567898.1 KL567913.1 KL568185.1 KL567892.1 KL567954.1 KL568053.1 KL568191.1 KL568322.1 KL568314.1 KL568076.1 KL567917.1 KL568092.1 KL567902.1 KL568420.1 KL567903.1 KL568122.1 KL567957.1 AABR07024192.1 KL568436.1 KL567906.1 KL568134.1 AABR07024041.1 KL568193.1 AABR07046231.1 KL568439.1 KL568054.1 KL568387.1 KL568085.1 KL567893.1 KL568457.1 KL568330.1 KL568481.1 KL568399.1 KL567894.1 KL568196.1 KL568130.1 KL568475.1 KL568006.1 KL568018.1 KL568133.1 KL568109.1 KL568307.1 KL568400.1 KL567882.1 KL568406.1 KL568010.1 KL567881.1 KL568068.1 KL568024.1 KL567937.1 KL568059.1 KL568072.1 KL568485.1 KL567933.1 KL567958.1 KL568009.1 KL568228.1 KL567895.1 KL568321.1 KL567975.1 KL567993.1 KL567965.1 AABR07024263.1 KL568016.1 KL568078.1 KL568391.1 KL568359.1 KL568295.1 KL567972.1 KL567924.1 KL567944.1 KL567986.1 KL568427.1 KL568507.1 KL568356.1 KL568100.1 KL567976.1 KL568404.1 KL568511.1 AABR07024032.1 KL568021.1 KL567889.1 KL568020.1 KL568011.1 KL568343.1 KL568389.1 KL568046.1 KL568088.1 KL568512.1 KL568490.1 KL568488.1 KL567911.1 KL568478.1 KL568052.1 KL568002.1 KL568041.1 KL568129.1 KL568012.1 KL568394.1 KL568376.1 KL568042.1 KL568248.1 KL567910.1 KL568403.1 KL568303.1 KL568245.1 KL568492.1 KL568169.1 KL568349.1 AABR07024044.1 KL567938.1 AABR07024071.1 KL568176.1 KL568411.1 KL568044.1 AABR07024206.1 KL568474.1 KL568312.1 KL567943.1 KL568087.1 AABR07024350.1 KL568270.1 KL567915.1 KL567884.1 KL568025.1 KL568007.1 KL568049.1 KL567916.1 KL567995.1 KL568219.1 AABR07024122.1 KL568036.1 KL567931.1 AABR07024286.1 KL567920.1 KL567999.1 KL568081.1 KL568111.1 KL567953.1 KL568039.1 KL568212.1 KL568089.1 KL568113.1 KL568272.1 KL568048.1 KL568038.1 KL568135.1 KL567880.1 KL568104.1 KL567970.1 KL568329.1 KL567945.1 KL568017.1 KL568066.1 KL568058.1 AABR07046136.1 KL567974.1 KL568317.1 KL568316.1 KL567886.1 KL568351.1 KL568486.1 AABR07024262.1 KL568084.1 KL568026.1 KL568380.1 KL568063.1 KL568015.1 KL568428.1 KL568075.1 KL568198.1 KL567921.1 KL568251.1 KL568441.1 KL567973.1 KL568215.1 KL567982.1 KL567935.1 KL567971.1 KL568121.1 KL568033.1 KL567968.1 KL568476.1 KL567936.1 KL567987.1 KL568218.1 KL568323.1 KL567980.1 KL567900.1 KL568390.1 KL568340.1 KL568132.1 KL568498.1 KL567967.1 KL568224.1 AABR07024394.1 KL568117.1 KL568331.1 KL568131.1 KL568437.1 KL567927.1 KL568517.1 KL567983.1 KL568499.1 KL568070.1 KL568073.1 KL568199.1 KL568384.1 KL568456.1 KL568170.1 KL568099.1 KL567952.1 AABR07024119.1 KL568285.1 KL568174.1 KL568396.1 KL568373.1 KL568093.1 KL568455.1 KL568484.1 KL568190.1 KL568514.1 KL568448.1 KL567904.1 KL568425.1 KL568030.1 KL568495.1 KL567963.1 KL568505.1 KL568283.1 KL568040.1 KL567909.1 KL567977.1 KL568326.1 KL567966.1 AABR07024264.1 KL567998.1 KL568167.1 KL568225.1 KL567899.1 KL568173.1 KL568341.1 KL568074.1 AABR07024040.1 KL568032.1 AABR07024150.1 KL567961.1 KL568304.1 KL568300.1 KL568335.1 KL567923.1 KL568127.1 KL568358.1 KL568112.1 KL568266.1 KL568091.1 KL568077.1 KL567912.1 AABR07024205.1 KL567896.1 KL568014.1 KL568028.1 KL567932.1 KL568110.1 KL568031.1 KL567942.1 KL568179.1 KL568383.1 KL568332.1 KL568095.1 KL568086.1 KL568022.1 KL568288.1 KL568180.1 KL568367.1 KL568223.1 KL568035.1 KL568235.1 KL568202.1 KL568055.1 KL568382.1 AABR07024321.1 KL568037.1 KL567885.1 KL567934.1 AABR07046230.1 KL567950.1 KL568493.1 KL568306.1 KL568208.1 KL567964.1 KL568310.1 KL568197.1 KL568226.1 KL567962.1 KL568013.1 KL568393.1 KL568263.1 KL568313.1 KL568071.1 KL567888.1 KL568062.1 KL568293.1 KL568195.1 KL568274.1 AABR07024292.1 KL568246.1 KL568177.1 KL568516.1 KL568128.1 KL567887.1 KL568234.1 KL568106.1 KL567922.1 KL568201.1 KL568297.1 KL568254.1 KL568268.1 KL567930.1 KL568271.1 KL568286.1 KL568518.1 KL568290.1 KL567918.1 KL568281.1 KL568333.1 KL568469.1 KL568289.1 KL568200.1 AABR07024101.1 KL568497.1 KL568513.1 KL568299.1 KL568256.1 KL568098.1 KL568043.1 KL567985.1 KL568280.1 KL568061.1 KL567949.1 KL568242.1 KL568090.1 AABR07024046.1 AABR07051059.1 KL568107.1 AABR07023990.1 KL568103.1 AABR07024118.1 KL567897.1 KL567960.1 KL568363.1 KL567879.1 AABR07024105.1 KL567928.1 AABR07024115.1 KL568259.1 KL567969.1 AABR07024067.1 AABR07024222.1 KL568120.1 KL568138.1 KL568372.1 AABR07045975.1 KL568047.1 AABR07051155.1 KL568123.1 KL568237.1 KL568260.1 KL568186.1 AABR07046142.1 AABR07024123.1 KL567978.1 KL568273.1 AABR07050784.1 KL568347.1 AABR07024193.1 KL568004.1 KL568365.1 KL568267.1 AABR07024228.1 AABR07024187.1 AABR07024267.1 AABR07022707.1 KL568247.1 KL568496.1 KL568339.1 KL568327.1 KL568392.1 KL567883.1 AABR07022511.1 KL568287.1 KL568000.1 AABR07050986.1 KL568444.1 KL568124.1 KL568309.1 KL567905.1 KL568060.1 AABR07024159.1 KL568094.1 AABR07050910.1 KL568459.1 KL568214.1 AABR07046212.1 AABR07046159.1 KL568192.1 KL568034.1 KL568233.1 KL567948.1 KL568311.1 KL568183.1 KL568401.1 KL568213.1 KL568069.1 KL568416.1 AABR07045689.1 AABR07024290.1 AABR07051090.1 KL568368.1 AABR07046158.1 KL568220.1 AABR07046187.1 KL568207.1 AABR07046222.1 KL568082.1 KL567955.1 AABR07055056.1 KL568056.1 KL568291.1 KL568083.1 KL567929.1 KL568050.1 AABR07024390.1 KL567996.1 AABR07022630.1 KL568205.1 AABR07022737.1 AABR07046433.1 KL568426.1 KL568503.1 AABR07046319.1 KL568227.1 KL568305.1 KL568301.1 KL568509.1 KL568320.1 KL568027.1 AABR07045776.1 KL568008.1 KL568172.1 KL568294.1 AABR07045688.1 KL568105.1 KL568232.1 AABR07050911.1 KL568194.1 AABR07024397.1 KL568378.1 AABR07024165.1 KL568258.1 AABR07045766.1 AABR07046296.1 KL568491.1 AABR07051126.1 KL568352.1 KL568188.1 AABR07050959.1 KL568515.1 AABR07022923.1 AABR07024399.1 AABR07023068.1 KL567992.1 KL568252.1 KL567990.1 AABR07050727.1 AABR07023168.1 KL568171.1 KL568231.1 KL568178.1 KL568435.1 KL567890.1 KL568221.1 KL568407.1 KL568065.1 AABR07022620.1 KL568217.1 AABR07024045.1 KL568262.1 KL568168.1 AABR07050989.1 KL568244.1 KL568241.1 AABR07022654.1 AABR07050887.1 KL568371.1 KL568182.1 KL568502.1 KL568175.1 AABR07022594.1 AABR07023325.1 AABR07024229.1 AABR07046267.1 KL568253.1 KL567925.1 KL568500.1 KL568240.1 KL568184.1 KL568278.1 KL568308.1 KL568243.1 KL568210.1 KL568222.1 KL568276.1 KL568429.1 KL568239.1 KL568385.1 KL568115.1 AABR07024323.1 KL568209.1 KL567989.1 KL568023.1 AABR07024195.1 KL568489.1 KL568302.1 KL568264.1 AABR07024204.1 AABR07024427.1 KL568261.1 KL568255.1 KL568284.1 KL568236.1 KL568203.1 AABR07024428.1 KL568211.1 KL568375.1 KL568229.1 KL568216.1 KL568501.1 AABR07046137.1 KL568189.1 KL568277.1 AABR07055065.1 KL568230.1 KL568338.1 KL568369.1 AABR07023968.1 KL568504.1 KL568386.1 AABR07024265.1 AABR07050825.1 KL568250.1 AABR07024293.1 AABR07024332.1 KL568249.1 AABR07046507.1 AABR07050815.1 KL568298.1 AABR07051099.1 KL568506.1 KL568494.1 KL568508.1 AABR07046249.1 KL568269.1 AABR07022926.1 KL568181.1 KL568265.1 KL568296.1 AABR07055042.1 AABR07024331.1 AABR07046563.1 AABR07022299.1 KL568422.1 AABR07023052.1 KL568238.1 AABR07024178.1 AABR07022871.1 AABR07024196.1 AABR07022830.1 AABR07024374.1 KL567991.1 KL568275.1 AABR07024221.1 AABR07022573.1 KL568279.1 AABR07046243.1 AABR07024269.1 AABR07046211.1 AABR07046150.1 KL568204.1 AABR07022857.1 AABR07024322.1 AABR07024080.1 AABR07024160.1 AABR07046014.1 AABR07022517.1 AABR07024182.1 AABR07050811.1 AABR07024125.1 AABR07024138.1 AABR07046186.1 AABR07022258.1 AABR07022378.1 AABR07024114.1 AABR07024103.1 AABR07024007.1 AABR07050916.1 AABR07024395.1 AABR07022405.1 AABR07022334.1 AABR07022616.1 AABR07023031.1 AABR07022626.1 AABR07022627.1 KL568206.1 AABR07022608.1 AABR07023445.1 AABR07022973.1 AABR07022342.1 AABR07023067.1 AABR07022531.1 AABR07045742.1 AABR07022905.1 AABR07050689.1 AABR07022404.1 AABR07024223.1 AABR07023024.1 AABR07022350.1 AABR07024383.1 AABR07022934.1 AABR07050912.1 AABR07050985.1 AABR07022759.1 AABR07022729.1 AABR07024396.1 AABR07022336.1 AABR07022359.1 AABR07024371.1 AABR07022264.1 AABR07022304.1 AABR07045681.1 AABR07024319.1 AABR07022672.1 AABR07022279.1 AABR07045743.1 AABR07022284.1 AABR07024227.1 AABR07024289.1 AABR07022416.1 AABR07022676.1 AABR07022287.1 AABR07046405.1 AABR07022684.1 AABR07024294.1 AABR07022265.1 AABR07024393.1 AABR07022907.1 AABR07024092.1 AABR07022259.1 KL567956.1 AABR07024349.1 AABR07024005.1 AABR07024043.1 AABR07022828.1 AABR07022519.1 AABR07022619.1 AABR07024189.1 AABR07046182.1 AABR07022760.1 AABR07024006.1 AABR07023321.1 KL568119.1 AABR07022335.1 AABR07022266.1 AABR07022921.1 AABR07023534.1 AABR07022518.1 AABR07022881.1 AABR07022993.1 AABR07023006.1 AABR07022367.1 AABR07024433.1 AABR07024351.1 AABR07022263.1 AABR07022990.1 AABR07022827.1 AABR07046248.1 AABR07022884.1 AABR07022629.1 AABR07024320.1 AABR07022559.1 AABR07022239.1 AABR07022428.1 AABR07024177.1 AABR07024266.1 AABR07023560.1 AABR07050894.1 AABR07022726.1 AABR07022738.1 AABR07024121.1 AABR07050661.1 AABR07022300.1 AABR07051089.1 AABR07024432.1 AABR07022880.1 AABR07022957.1 AABR07051098.1 KL567979.1 AABR07024304.1 AABR07022922.1 AABR07046283.1 AABR07022337.1 AABR07022671.1 AABR07022257.1 AABR07022783.1 AABR07022425.1 AABR07024194.1 AABR07022273.1 AABR07022504.1 AABR07022368.1 AABR07024226.1 AABR07022829.1 AABR07050850.1 AABR07022621.1 AABR07022628.1 AABR07022431.1 AABR07022478.1 AABR07046034.1 AABR07022753.1 AABR07023014.1 AABR07022994.1 AABR07022645.1 AABR07046019.1 AABR07022288.1 AABR07022802.1 AABR07051154.1 AABR07023967.1 AABR07022920.1 AABR07023949.1 AABR07045875.1 AABR07024100.1 AABR07022479.1 AABR07022532.1 AABR07023324.1 AABR07022793.1 AABR07022301.1 AABR07022906.1 AABR07022514.1 AABR07022298.1 AABR07022333.1 AABR07022456.1 AABR07022885.1 AABR07022944.1 AABR07024268.1"
                                                           ]})

# Define general variables
genome_used = (str(config["genome"])).lower()
blacklist = (genomes_table[genomes_table['genome_id']==genome_used]).blacklist.iloc[0]
genomeSize = (genomes_table[genomes_table['genome_id']==genome_used]).effective_genomeSize.iloc[0]
ignore_for_normalization = (genomes_table[genomes_table['genome_id']==genome_used]).ignore_for_normalization.iloc[0]

if ((eval(str(config["paired_end"])) == True)):
    read_extension = "--extendReads"
else:
    read_extension = "--extendReads "+str(config["fragment_length"])


if ((eval(str(config["use_macs3"])) == True)):
    macs_version = "macs3"
else:
    macs_version = "macs2"


### working directory
home_dir = os.path.join(config["output_directory"],"")
shell('mkdir -p {home_dir}')
workdir: home_dir


### get the unique samples names and other variables
# loading the sample table
sample_metadata = pd.read_csv(str(config["sample_config_table"]),  sep='\t+', engine='python')   # target_id | input_id | broad
sample_metadata = sample_metadata.iloc[:,0:3].set_axis(['target_id', 'input_id', 'broad'], axis=1, inplace=False)
TARGETNAMES = list(numpy.unique(list(sample_metadata.target_id)))
INPUTNAMES = list(numpy.unique(list(sample_metadata.input_id)))
SAMPLENAMES = list(numpy.unique(TARGETNAMES + INPUTNAMES))


# Get bam list
if not (os.path.exists(config["runs_directory"])):
    os.system("printf '\033[1;31m\\n!!! *runs_directory* does not exist !!!\\n\\n\033[0m'")
else:
    BAMS = next(os.walk(config["runs_directory"]))[2]
    RUNNAMES = numpy.unique([re.sub(rf"{config['bam_suffix']}$", "", i) for i in BAMS])



### MACS2 or MACS3?
if ((eval(str(config["use_macs3"])) == True)):
    PEAKCALLER = "macs3"
else:
    PEAKCALLER = "macs2"

### other generic variable
if ((eval(str(config["remove_duplicates"])) == True)):
    DUP = "dedup"
else:
    DUP = "mdup"



### Optional analysis outputs
if ((eval(str(config["perform_plotFingerprint"])) == True)):
    plotFingerprint_results = expand("05_Quality_controls_and_statistics/plotFingerprint/{target}_fingerPrinting_plot.pdf", target = TARGETNAMES)
else:
    plotFingerprint_results = []

if ((eval(str(config["perform_fragmentSizeDistribution"])) == True) & (eval(str(config["paired_end"])) == True)):
    fragmentSizeDistribution_results = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_metrics.txt"
else:
    fragmentSizeDistribution_results = []

if (len(SAMPLENAMES) > 2):
    PCA_wholeGenome_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.1-2_heatmap_wholeGenome.pdf"
    PCA_wholeGenome_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.2-3_heatmap_wholeGenome.pdf"
else:
    PCA_wholeGenome_12 = []
    PCA_wholeGenome_23 = []


if (len(TARGETNAMES) > 2):
    PCA_atPeaks_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.1-2_heatmap_atPeaks.pdf"
    PCA_atPeaks_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.2-3_heatmap_atPeaks.pdf"
else:
    PCA_atPeaks_12 = []
    PCA_atPeaks_23 = []


if (len(SAMPLENAMES) > 1):
    correlation_heatmap_wholeGenome_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome.pdf"
    correlation_heatmap_wholeGenome_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome.pdf"
else:
    correlation_heatmap_wholeGenome_pearson = []
    correlation_heatmap_wholeGenome_spearman = []


if (len(TARGETNAMES) > 1):
    correlation_heatmap_atPeaks_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_pearson.correlation_heatmap_atPeaks.pdf"
    correlation_heatmap_atPeaks_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_spearman.correlation_heatmap_atPeaks.pdf"
else:
    correlation_heatmap_atPeaks_pearson = []
    correlation_heatmap_atPeaks_spearman = []


if ((eval(str(config["paired_end"])) == True)):
    peaks = expand("04_Called_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES)
else:
    peaks = expand("04_Called_peaks/{target}.filtered.BAM_peaks.xls", target = TARGETNAMES)


### Generation of global wildcard_constraints
# Function to handle the values for the wilcards
def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each
    value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    """
    if isinstance(values, str):
            raise ValueError("constraint_to(): Expected a list, got str instead")
    return "({})".format("|".join(values))

wildcard_constraints:
    SAMPLE = constraint_to(SAMPLENAMES),
    TARGET = constraint_to(TARGETNAMES),
    INPUT = constraint_to(INPUTNAMES)


ruleorder: fastQC_filtered_BAM > normalized_bigWig > raw_bigWig

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================
# Function to run all funtions
if (set(SAMPLENAMES) <= set(RUNNAMES)):
    rule AAA_initialization:
        input:
            fastqc_bam_zip = expand(os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"])), sample = SAMPLENAMES),
            plotFingerprint_results = plotFingerprint_results,
            fragmentSizeDistribution_results = fragmentSizeDistribution_results,
            normalized_bigWig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), sample = SAMPLENAMES),
            raw_bigWig = expand(os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_raw.coverage_bs", str(config["bigWig_binSize"]), ".bw"])), sample = SAMPLENAMES),
            correlation_heatmap_wholeGenome_pearson = correlation_heatmap_wholeGenome_pearson,
            correlation_heatmap_wholeGenome_spearman = correlation_heatmap_wholeGenome_spearman,
            PCA_wholeGenome_12 = PCA_wholeGenome_12,
            PCA_wholeGenome_23 = PCA_wholeGenome_23,
            peaks = peaks,
            multiqc_report = "05_Quality_controls_and_statistics/multiQC/multiQC_report.html",
            correlation_heatmap_atPeaks_pearson = correlation_heatmap_atPeaks_pearson,
            correlation_heatmap_atPeaks_spearman = correlation_heatmap_atPeaks_spearman,
            PCA_atPeaks_12 = PCA_atPeaks_12,
            PCA_atPeaks_23 = PCA_atPeaks_23,
            aggregated_qc = "05_Quality_controls_and_statistics/peaks_stats/all_samples_FRiP_report.tsv"
        shell:
            """
            printf '\033[1;36mPipeline ended!\\n\033[0m'
            """
else:
    missing_samples = '\\n  - '.join(list(set(SAMPLENAMES) - set(RUNNAMES)))
    os.system("printf '\033[1;31m\\n!!! Not all bam files are avalable in the input directory !!!\\n\\nPlease provide files for:\\n  - "+missing_samples+"\\n\\n\033[0m'")

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================

if (eval(str(config["skip_bam_filtering"])) == False):
    rule MAPQ_filter:
        input:
            source_bam = os.path.join(config["runs_directory"], ''.join(["{SAMPLE}", config["bam_suffix"]]))
        output:
            bam_mapq_only = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), ".bam"]))),
            bam_mapq_only_sorted = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam"]))),
            bam_mapq_only_sorted_index = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam.bai"])))
        params:
            sample = "{SAMPLE}",
            MAPQ_threshold = config["MAPQ_threshold"]
        threads:
            math.ceil(workflow.cores/len(SAMPLENAMES))
        shell:
            """
            printf '\033[1;36m{params.sample}: filtering MAPQ and re-indexing...\\n\033[0m'

            samtools view -@ {threads} -h -q {params.MAPQ_threshold} {input.source_bam} -o {output.bam_mapq_only}

            samtools sort -@ {threads} {output.bam_mapq_only} -o {output.bam_mapq_only_sorted}
            samtools index -@ {threads} -b {output.bam_mapq_only_sorted} {output.bam_mapq_only_sorted_index}
            """


    if ((eval(str(config["paired_end"])) == True) & (eval(str(config["umi_present"])) == True)):
        rule gatk4_markdups_umiAware:
            input:
                bam_mapq_only_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam"])),
                bam_mapq_only_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam.bai"]))
            output:
                bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
                bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
                umi_metrics = "01_BAM_filtered/umi_metrics/{SAMPLE}_UMI_metrics.txt",
                dup_metrics = "01_BAM_filtered/MarkDuplicates_metrics/{SAMPLE}_MarkDuplicates_metrics.txt",
                flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"]))
            params:
                remove_duplicates = (str(config["remove_duplicates"])).lower(),
                sample = "{SAMPLE}"
            log:
                out = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.out",
                err = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.err"
            threads:
                workflow.cores
            shell:
                """
                printf '\033[1;36m{params.sample}: UMI-aware gatk MarkDuplicates...\\n\033[0m'

                mkdir -p 01_BAM_filtered/umi_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_logs
                mkdir -p 01_BAM_filtered/flagstat

                gatk UmiAwareMarkDuplicatesWithMateCigar \
                --INPUT {input.bam_mapq_only_sorted} \
                --OUTPUT {output.bam_mdup} \
                --REMOVE_DUPLICATES {params.remove_duplicates} \
                --MAX_EDIT_DISTANCE_TO_JOIN 1 \
                --UMI_METRICS_FILE {output.umi_metrics} \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --UMI_TAG_NAME RX \
                --CREATE_INDEX true \
                --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

                samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
                """
    else: # Single-end/no-UMI dedup
        rule gatk4_markdups:
            input:
                bam_mapq_only_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam"])),
                bam_mapq_only_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam.bai"]))
            output:
                bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
                bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
                dup_metrics = "01_BAM_filtered/MarkDuplicates_metrics/{SAMPLE}_MarkDuplicates_metrics.txt",
                flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"]))
            params:
                remove_duplicates = (str(config["remove_duplicates"])).lower(),
                sample = "{SAMPLE}"
            log:
                out = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.out",
                err = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.err"
            threads:
                workflow.cores
            shell:
                """
                printf '\033[1;36m{params.sample}: 'standard' gatk MarkDuplicates...\\n\033[0m'

                mkdir -p 01_BAM_filtered/MarkDuplicates_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_logs
                mkdir -p 01_BAM_filtered/flagstat

                gatk MarkDuplicatesWithMateCigar \
                --INPUT {input.bam_mapq_only_sorted} \
                --OUTPUT {output.bam_mdup} \
                --REMOVE_DUPLICATES {params.remove_duplicates} \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --CREATE_INDEX true \
                --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

                samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
                """
else:
    rule bam_link__skip_filtering:
        input:
            source_bam = os.path.join(config["runs_directory"], ''.join(["{SAMPLE}", config["bam_suffix"]]))
        output:
            bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"]))
        params:
            sample = "{SAMPLE}"
        threads:
            math.ceil(workflow.cores/len(SAMPLENAMES))
        shell:
            """
            printf '\033[1;36m{params.sample} (skip filtering): linking bam, indexing and computing flagstat...\\n\033[0m'

            mkdir -p 01_BAM_filtered/flagstat

            BAM_REAL=$(realpath {input.source_bam})
            ln -s $BAM_REAL {output.bam_mdup}
            samtools index -@ {threads} -b {output.bam_mdup} {output.bai_mdup}

            samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
            """



rule fastQC_filtered_BAM:
    input:
        bam_mapq = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        bai_mapq = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        html = os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.html"])),
        zip = os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"]))
    params:
        fastQC_BAMs_outdir = os.path.join("02_fastQC_on_BAM_filtered/"),
        sample = "{SAMPLE}"
    threads:
        math.ceil(workflow.cores/len(SAMPLENAMES))
    shell:
        """
        mkdir -p 02_fastQC_on_BAM_filtered

        printf '\033[1;36m{params.sample}: Performing fastQC on deduplicated bam...\\n\033[0m'
        fastqc -t {threads} --outdir {params.fastQC_BAMs_outdir} {input.bam_mapq}
        """

# ------------------------------------------------------------------------------

rule plotFingerprint:
    input:
        target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        lorenz_curve_pdf = os.path.join("05_Quality_controls_and_statistics/plotFingerprint/{TARGET}_fingerPrinting_plot.pdf"),
        quality_metrics = os.path.join("05_Quality_controls_and_statistics/plotFingerprint/quality_metrics/{TARGET}_fingerPrinting_quality_metrics.txt")
    params:
        sample = "{TARGET}",
        sample_config_table = config["sample_config_table"],
        input_suffix = "_mapq"+str(config["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
        read_extension = read_extension,
        blacklist = blacklist
    threads:
        math.ceil(workflow.cores/len(TARGETNAMES))
    log:
        out = "05_Quality_controls_and_statistics/plotFingerprint/logs/{TARGET}_fingerPrinting_log.out",
        err = "05_Quality_controls_and_statistics/plotFingerprint/logs/{TARGET}_fingerPrinting_log.err"
    shell:
        """
        printf '\033[1;36m{params.sample}: plotting fingerprint...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/plotFingerprint/logs
        mkdir -p 05_Quality_controls_and_statistics/plotFingerprint/quality_metrics

        INPUT_ID=$(grep {params.sample} {params.sample_config_table} | cut -f 2)

        plotFingerprint \
        -b {input.target_bam} \
        01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
        --JSDsample 01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
        -plot {output.lorenz_curve_pdf} \
        {params.read_extension} \
        --ignoreDuplicates \
        --outQualityMetrics {output.quality_metrics} \
        --labels {params.sample} ${{INPUT_ID}} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule fragmentSizeDistribution:
    input:
        all_bams = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), sample = SAMPLENAMES),
        all_bais = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])), sample = SAMPLENAMES)
    output:
        fragment_distribution_plot = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_plot.pdf",
        fragmentSize_metrics = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_metrics.txt",
        fragmentSize_RawFragmentLengths = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_RawFragmentLengths.tab"
    params:
        labels = ' '.join(SAMPLENAMES),
        blacklist = blacklist
    threads:
        workflow.cores
    log:
        out = "05_Quality_controls_and_statistics/fragmentSize_distribution/logs/fragmentSize_distribution_log.out",
        err = "05_Quality_controls_and_statistics/fragmentSize_distribution/logs/fragmentSize_distribution_log.err"
    shell:
        """
        printf '\033[1;36mPlotting fragment size distribution...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/fragmentSize_distribution/logs

        bamPEFragmentSize \
        --bamfiles {input.all_bams} \
        --binSize 1000000 \
        --blackListFileName {params.blacklist} \
        --samplesLabel {params.labels} \
        --histogram {output.fragment_distribution_plot} \
        --table {output.fragmentSize_metrics} \
        --outRawFragmentLengths {output.fragmentSize_RawFragmentLengths} \
        -p {threads} > {log.out} 2> {log.err}
        """

# ------------------------------------------------------------------------------

rule normalized_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        normalized_bigWig = os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])),
    params:
        sample = "{SAMPLE}",
        blacklist = blacklist,
        genomeSize = genomeSize,
        ignore_for_normalization = ignore_for_normalization,
        read_extension = read_extension,
        bw_binSize = config["bigWig_binSize"]
    threads:
        math.ceil(workflow.cores/len(SAMPLENAMES))
    log:
        out = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLE}_fragmentSize_distribution_log.out",
        err = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLE}_fragmentSize_distribution_log.err"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating RPGC normalized bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/RPGC_normalized/logs

        bamCoverage \
        -b {input.bam} \
        -o {output.normalized_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing RPGC \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        --ignoreDuplicates \
        {params.read_extension} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule raw_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        raw_bigWig = os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_raw.coverage_bs", str(config["bigWig_binSize"]), ".bw"])),
    params:
        sample = "{SAMPLE}",
        blacklist = blacklist,
        genomeSize = genomeSize,
        ignore_for_normalization = ignore_for_normalization,
        read_extension = read_extension,
        bw_binSize = config["bigWig_binSize"]
    threads:
        math.ceil(workflow.cores/len(SAMPLENAMES))
    log:
        out = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLE}_fragmentSize_distribution_log.out",
        err = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLE}_fragmentSize_distribution_log.err"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating raw coverage bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/raw_coverage/logs

        bamCoverage \
        -b {input.bam} \
        -o {output.raw_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing None \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        --ignoreDuplicates \
        {params.read_extension} \
        -p {threads} > {log.out} 2> {log.err}
        """

# ------------------------------------------------------------------------------

rule multiBigwigSummary_wholeGenome:
    input:
        all_norm_bigwig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), sample = SAMPLENAMES),
    output:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome.npz"
    params:
        labels = ' '.join(SAMPLENAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads:
        workflow.cores
    log:
        out = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/multiBigWigSummary_matrix_wholeGenome_log.out",
        err = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/multiBigWigSummary_matrix_wholeGenome_log.err"
    shell:
        """
        printf '\033[1;36mComputing multiBigwigSummary matrix (whole genome)...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs

        multiBigwigSummary bins \
        -b {input.all_norm_bigwig} \
        -o {output.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --binSize 1000 \
        --chromosomesToSkip {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule correlations_wholeGenome:
    input:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome.npz"
    output:
        correlation_heatmap_wholeGenome_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome.pdf",
        correlation_heatmap_wholeGenome_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome.pdf"
    params:
        labels = ' '.join(SAMPLENAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_pearson.correlation_heatmap_wholeGenome_log.out",
        err_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_pearson.correlation_heatmap_wholeGenome_log.err",
        out_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_spearman.correlation_heatmap_wholeGenome_log.out",
        err_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_spearman.correlation_heatmap_wholeGenome_log.err"
    shell:
        """
        printf '\033[1;36mPlotting sample correlations (whole genome)...\\n\033[0m'

        plotCorrelation \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --plotTitle 'Pearson correlation whole genome RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_wholeGenome_pearson} \
        --colorMap {params.heatmap_color} > {log.out_pearson} 2> {log.err_pearson}


        plotCorrelation \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --corMethod spearman \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --plotTitle 'Spearman correlation whole genome RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_wholeGenome_spearman} \
        --colorMap {params.heatmap_color} > {log.out_spearman} 2> {log.err_spearman}
        """



rule PCA_wholeGenome:
    input:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome.npz"
    output:
        PCA_wholeGenome_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.1-2_heatmap_wholeGenome.pdf",
        PCA_wholeGenome_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.2-3_heatmap_wholeGenome.pdf"
    params:
        labels = ' '.join(SAMPLENAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.1-2_heatmap_wholeGenome_log.out",
        err_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.1-2_heatmap_wholeGenome_log.err",
        out_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.2-3_heatmap_wholeGenome_log.out",
        err_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.2-3_heatmap_wholeGenome_log.err"
    shell:
        """
        printf '\033[1;36mPlotting PCA (whole genome)...\\n\033[0m'

        plotPCA \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --PCs 1 2 \
        --plotTitle 'PCA whole genome: PC1 vs PC2 (RPGC normalized coverage)' \
        --plotFile {output.PCA_wholeGenome_12} > {log.out_12} 2> {log.err_12}

        plotPCA \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --PCs 2 3 \
        --plotTitle 'PCA whole genome: PC2 vs PC3 (RPGC normalized coverage)' \
        --plotFile {output.PCA_wholeGenome_23} > {log.out_23} 2> {log.err_23}
        """

# ------------------------------------------------------------------------------

if ((eval(str(config["paired_end"])) == True)):
    rule macs_callpeak_PE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            input_bam_all = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input = INPUTNAMES),
        output:
            peaksPE = "04_Called_peaks/{TARGET}.filtered.BAMPE_peaks.xls"
        params:
            sample = "{TARGET}",
            macs_version = macs_version,
            sample_config_table = config["sample_config_table"],
            input_suffix = "_mapq"+str(config["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
            read_extension = read_extension,
            genomeSize = genomeSize,
            macs_qValue_cutoff = config["macs_qValue_cutoff"],
            blacklist = blacklist
        threads:
            math.ceil(workflow.cores/len(TARGETNAMES))
        log:
            out = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAMPE_log.out",
            err = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAMPE_log.err"
        shell:
            """
            printf '\033[1;36m{params.sample}: calling peaks ({params.macs_version})...\\n\033[0m'

            mkdir -p 04_Called_peaks/logs

            INPUT_ID=$(grep {params.sample} {params.sample_config_table} | cut -f 2)
            CALL_BROAD=$(grep {params.sample} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                BROAD=""
            else
                BROAD="--broad"
            fi

            {params.macs_version} callpeak \
            -t {input.target_bam} \
            -c 01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
            -f BAMPE \
            -g {params.genomeSize} \
            -q {params.macs_qValue_cutoff} \
            --keep-dup all \
            --outdir 04_Called_peaks \
            --name {params.sample}.filtered.BAMPE ${{BROAD}} > {log.err} 2> {log.out}
            """
else:
    rule phantom_SE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            input_bam_all = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input = INPUTNAMES)
        output:
            phantom = '04_Called_peaks/phantom/{TARGET}.phantom.spp.out'
        log:
            out = '04_Called_peaks/phantom/logs/{TARGET}.phantom.log'
        params:
            sample = "{TARGET}",
            sample_config_table = config["sample_config_table"],
            input_suffix = "_mapq"+str(config["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
            genomeSize = genomeSize,
            macs_qValue_cutoff = config["macs_qValue_cutoff"],
            blacklist = blacklist
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36m{params.sample}: calculating phantom peak...\\n\033[0m'
            mkdir -p 04_Called_peaks/phantom/logs

            INPUT_ID=$(grep {params.sample} {params.sample_config_table} | cut -f 2)

            ${{CONDA_PREFIX}}/bin/Rscript ${{CONDA_PREFIX}}/bin/run_spp.R -rf -c='{input.target_bam}' -i="01_BAM_filtered/${{INPUT_ID}}{params.input_suffix}" -savp -out={output.phantom} &> {log.out}
            """



    rule fragment_length:
        input:
            phantom = '04_Called_peaks/phantom/{TARGET}.phantom.spp.out'
        output:
            fragment_length_phanthom = temp('04_Called_peaks/phantom/{TARGET}.fragment_length')
        shell:
            """
            awk '{{print $3}}' < {input.phantom} | tr ',' '\\t' | awk '{{if($1!=0) print $1; else print $2}}' > {output.fragment_length_phanthom}
            """



    rule macs_callpeak_SE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            input_bam_all = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input = INPUTNAMES),
            phantom = '04_Called_peaks/phantom/{TARGET}.fragment_length'
        output:
            peaksSE = "04_Called_peaks/{TARGET}.filtered.BAM_peaks.xls"
        params:
            sample = "{TARGET}",
            macs_version = macs_version,
            sample_config_table = config["sample_config_table"],
            input_suffix = "_mapq"+str(config["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
            genomeSize = genomeSize,
            macs_qValue_cutoff = config["macs_qValue_cutoff"],
            blacklist = blacklist
        threads:
            math.ceil(workflow.cores/len(TARGETNAMES))
        log:
            out = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAM_log.out",
            err = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAM_log.err"
        shell:
            """
            printf '\033[1;36m{params.sample}: calling peaks ({params.macs_version})...\\n\033[0m'

            mkdir -p 04_Called_peaks/logs

            INPUT_ID=$(grep {params.sample} {params.sample_config_table} | cut -f 2)
            CALL_BROAD=$(grep {params.sample} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                BROAD=""
            else
                BROAD="--broad"
            fi

            {params.macs_version} callpeak \
            -t {input.target_bam} \
            -c 01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
            -f BAM \
            -g {params.genomeSize} \
            --nomodel \
            -q {params.macs_qValue_cutoff} \
            --outdir 04_Called_peaks \
            --name {params.sample}.filtered.BAM \
            --extsize $(cat {input.phantom}) ${{BROAD}} > {log.err} 2> {log.out}
            """

# ------------------------------------------------------------------------------

if (eval(str(config["skip_bam_filtering"])) == False):
    picard_metrics_file = expand("01_BAM_filtered/MarkDuplicates_metrics/{sample}_MarkDuplicates_metrics.txt", sample = SAMPLENAMES)
    picard_metrics_dir = "01_BAM_filtered/MarkDuplicates_metrics"
else:
    picard_metrics_file = []
    picard_metrics_dir = []


if ((eval(str(config["paired_end"])) == True)):
    rule multiQC_PE:
        input:
            fastqc = expand(os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"])), sample = SAMPLENAMES),
            picard_metrics = picard_metrics_file,
            flagstat = expand(os.path.join("01_BAM_filtered/flagstat/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"])), sample = SAMPLENAMES),
            peaks = expand("04_Called_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES)
        output:
            multiqc_report = "05_Quality_controls_and_statistics/multiQC/multiQC_report.html"
        params:
            out_directory = "05_Quality_controls_and_statistics/multiQC/",
            multiqc_report_name = "multiQC_report.html",
            picard_metrics_dir = picard_metrics_dir
        threads: 1
        log:
            out = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.out",
            err = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.err"
        shell:
            """
            mkdir -p 05_Quality_controls_and_statistics/multiQC/

            $CONDA_PREFIX/bin/multiqc -f \
            -o {params.out_directory} \
            -n {params.multiqc_report_name} \
            --dirs 02_fastQC_on_BAM_filtered {params.picard_metrics_dir} 01_BAM_filtered/flagstat 04_Called_peaks > {log.err} 2> {log.out}
            """
else:
    rule multiQC_SE:
        input:
            fastqc = expand(os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"])), sample = SAMPLENAMES),
            picard_metrics = picard_metrics_file,
            flagstat = expand(os.path.join("01_BAM_filtered/flagstat/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"])), sample = SAMPLENAMES),
            peaks = expand("04_Called_peaks/{target}.filtered.BAM_peaks.xls", target = TARGETNAMES),
            phanthom = expand('04_Called_peaks/phantom/{target}.phantom.spp.out', target = TARGETNAMES)
        output:
            multiqc_report = "05_Quality_controls_and_statistics/multiQC/multiQC_report.html"
        params:
            out_directory = "05_Quality_controls_and_statistics/multiQC/",
            multiqc_report_name = "multiQC_report.html",
            picard_metrics_dir = picard_metrics_dir
        threads: 1
        log:
            out = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.out",
            err = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.err"
        shell:
            """
            printf '\033[1;36mGenerating multiQC report...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/multiQC/

            $CONDA_PREFIX/bin/multiqc -f \
            -o {params.out_directory} \
            -n {params.multiqc_report_name} \
            --dirs 02_fastQC_on_BAM_filtered {params.picard_metrics_dir} 01_BAM_filtered/flagstat 04_Called_peaks > {log.err} 2> {log.out}
            """

# ------------------------------------------------------------------------------

if ((eval(str(config["paired_end"])) == True)):
    rule merge_all_peaks_PE:
        input:
            peaks = expand("04_Called_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES)
        output:
            concat_peaks = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat.bed"),
            concat_peaks_sorted = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat_sorted.bed"),
            merged_peaks_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
        params:
            labels = ' '.join(SAMPLENAMES),
            blacklist = blacklist,
            ignore_for_normalization = ignore_for_normalization
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36mMerging all peaks...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_atPeaks/

            PEAK_LIST=$(ls 04_Called_peaks/*Peak | grep -v gapped)

            for i in ${{PEAK_LIST}}
            do
                cut -f 1-6 $i >> {output.concat_peaks}
            done

            sort -V -k1,1 -k2,2 {output.concat_peaks} > {output.concat_peaks_sorted}

            bedtools merge -i {output.concat_peaks_sorted} | sort -V -k1,1 -k2,2 > {output.merged_peaks_sorted}
            """
else:
    rule merge_all_peaks_SE:
        input:
            peaks = expand("04_Called_peaks/{target}.filtered.BAM_peaks.xls", target = TARGETNAMES)
        output:
            concat_peaks = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat.bed"),
            concat_peaks_sorted = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat_sorted.bed"),
            merged_peaks_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
        params:
            labels = ' '.join(SAMPLENAMES),
            blacklist = blacklist,
            ignore_for_normalization = ignore_for_normalization
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36mMerging all peaks...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_atPeaks/

            PEAK_LIST=$(ls 04_Called_peaks/*Peak | grep -v gapped)

            for i in ${{PEAK_LIST}}
            do
                cut -f 1-6 $i >> {output.concat_peaks}
            done

            sort -V -k1,1 -k2,2 {output.concat_peaks} > {output.concat_peaks_sorted}

            bedtools merge -i {output.concat_peaks_sorted} | sort -V -k1,1 -k2,2 > {output.merged_peaks_sorted}
            """



rule multiBigwigSummary_atPeaks:
    input:
        all_norm_bigwig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{target}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), target = TARGETNAMES),
        merged_peaks_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
    output:
        multiBigWig_matrix_atPeaks = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.npz"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads:
        workflow.cores
    log:
        out = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/multiBigWigSummary_matrix_atPeaks_log.out",
        err = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/multiBigWigSummary_matrix_atPeaks_log.err"
    shell:
        """
        printf '\033[1;36mComputing multiBigwigSummary matrix (at peaks)...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs

        multiBigwigSummary BED-file \
        --BED {input.merged_peaks_sorted} \
        -b {input.all_norm_bigwig} \
        -o {output.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --binSize 1000 \
        --chromosomesToSkip {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule correlations_atPeaks:
    input:
        multiBigWig_matrix_atPeaks = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.npz"
    output:
        correlation_heatmap_atPeaks_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_pearson.correlation_heatmap_atPeaks.pdf",
        correlation_heatmap_atPeaks_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_spearman.correlation_heatmap_atPeaks.pdf"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_pearson.correlation_heatmap_atPeaks_log.out",
        err_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_pearson.correlation_heatmap_atPeaks_log.err",
        out_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_spearman.correlation_heatmap_atPeaks_log.out",
        err_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_spearman.correlation_heatmap_atPeaks_log.err"
    shell:
        """
        printf '\033[1;36mPlotting sample correlations (at peaks)...\\n\033[0m'

        plotCorrelation \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --plotTitle 'Pearson correlation at peaks RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_atPeaks_pearson} \
        --colorMap {params.heatmap_color} > {log.out_pearson} 2> {log.err_pearson}


        plotCorrelation \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --corMethod spearman \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --plotTitle 'Spearman correlation at peaks RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_atPeaks_spearman} \
        --colorMap {params.heatmap_color} > {log.out_spearman} 2> {log.err_spearman}
        """



rule PCA_atPeaks:
    input:
        multiBigWig_matrix_atPeaks = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.npz"
    output:
        PCA_atPeaks_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.1-2_heatmap_atPeaks.pdf",
        PCA_atPeaks_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.2-3_heatmap_atPeaks.pdf"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads: 1
    log:
        out_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.1-2_heatmap_atPeaks_log.out",
        err_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.1-2_heatmap_atPeaks_log.err",
        out_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.2-3_heatmap_atPeaks_log.out",
        err_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.2-3_heatmap_atPeaks_log.err"
    shell:
        """
        printf '\033[1;36mPlotting PCA (at peaks)...\\n\033[0m'

        plotPCA \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --PCs 1 2 \
        --plotTitle 'PCA at peaks: PC1 vs PC2 (RPGC normalized coverage)' \
        --plotFile {output.PCA_atPeaks_12} > {log.out_12} 2> {log.err_12}

        plotPCA \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --PCs 2 3 \
        --plotTitle 'PCA at peaks: PC2 vs PC3 (RPGC normalized coverage)' \
        --plotFile {output.PCA_atPeaks_23} > {log.out_23} 2> {log.err_23}
        """

# ------------------------------------------------------------------------------

if ((eval(str(config["paired_end"])) == True)):
    rule MACS_peak_QC_PE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
	        peaks = "04_Called_peaks/{TARGET}.filtered.BAMPE_peaks.xls"
        output:
	        qc = temp("05_Quality_controls_and_statistics/peaks_stats/{TARGET}.filtered.BAMPE_peaks.qc.txt")
        params:
            sample_config_table = config["sample_config_table"],
            peak_prefix = "04_Called_peaks/{TARGET}.filtered.BAMPE_peaks",
            blacklist = blacklist,
            target = "{TARGET}",
            genomeSize = genomeSize
        threads:
            math.ceil(workflow.cores/len(TARGETNAMES))
        shell:
            """
            printf '\033[1;36m{params.target}: computing peak stats...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/peaks_stats/

            # define peak file
            CALL_BROAD=$(grep {params.target} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                CALLING_MODE="narrow"
                PEAK="{params.peak_prefix}.narrowPeak"
                PEAK_CHR="{params.peak_prefix}_chr.narrowPeak"
            else
                CALLING_MODE="broad"
                PEAK="{params.peak_prefix}.broadPeak"
                PEAK_CHR="{params.peak_prefix}_chr.broadPeak"
            fi

            # get the number of peaks
            peak_count=$(wc -l < $PEAK)

            # get the number of mapped reads
            mapped_reads=$(samtools view -c -F 4 {input.target_bam})

            # calculate the number of alignments overlapping the peaks
            # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
            reads_in_peaks=$(samtools view -@ {threads} -c -F 4 -L $PEAK {input.target_bam})

            # calculate Fraction of Reads In Peaks
            frip=$(bc -l <<< "$reads_in_peaks/$mapped_reads")

            # compute peak genome coverage
            peak_len=$(awk '{{total+=$3-$2}}END{{print total}}' $PEAK)
            genome_size={params.genomeSize}
            genomecov=$(bc -l <<< "$peak_len/$genome_size")

            # rounding fractions
            genomecov_round=$(printf "%.5f\n" "$genomecov")
            frip_round=$(printf "%.3f\n" "$frip")

            # write peak-based QC metrics to output file
            printf '{params.target}\\t'$CALLING_MODE'\\t'$peak_count'\\t'$frip_round'\\t'$genomecov_round'\\n' > {output.qc}

	        # add chr to peak files
            bedtools subtract -nonamecheck -a $PEAK -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > $PEAK_CHR
	        """

    rule aggregate_FRiP_PE:
        input:
            qc = expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAMPE_peaks.qc.txt", target = TARGETNAMES)
        output:
            aggregated_qc = "05_Quality_controls_and_statistics/peaks_stats/all_samples_FRiP_report.tsv"
        params:
            all_qc = ' '.join(expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAMPE_peaks.qc.txt", target = TARGETNAMES))
        threads: 1
        shell:
            """
            # print header of FRiP report
            printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
            cat {params.all_qc} >> {output.aggregated_qc}
            """

#*******************************************************************
else:
    rule MACS_peak_QC_SE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
	        peaks = "04_Called_peaks/{TARGET}.filtered.BAM_peaks.xls"
        output:
	        qc = temp("05_Quality_controls_and_statistics/peaks_stats/{TARGET}.filtered.BAM_peaks.qc.txt")
        params:
            sample_config_table = config["sample_config_table"],
            peak_prefix = "04_Called_peaks/{TARGET}.filtered.BAM_peaks",
            blacklist = blacklist,
            target = "{TARGET}",
            genomeSize = genomeSize
        threads:
            math.ceil(workflow.cores/len(TARGETNAMES))
        shell:
            """
            printf '\033[1;36m{params.target}: computing peak stats...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/peaks_stats/

            # define peak file
            CALL_BROAD=$(grep {params.target} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                CALLING_MODE="narrow"
                PEAK="{params.peak_prefix}.narrowPeak"
                PEAK_CHR="{params.peak_prefix}_chr.narrowPeak"
            else
                CALLING_MODE="broad"
                PEAK="{params.peak_prefix}.broadPeak"
                PEAK_CHR="{params.peak_prefix}_chr.broadPeak"
            fi

            # get the number of peaks
            peak_count=$(wc -l < $PEAK)

            # get the number of mapped reads
            mapped_reads=$(samtools view -c -F 4 {input.target_bam})

            # calculate the number of alignments overlapping the peaks
            # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
            reads_in_peaks=$(samtools view -@ {threads} -c -F 4 -L $PEAK {input.target_bam})

            # calculate Fraction of Reads In Peaks
            frip=$(bc -l <<< "$reads_in_peaks/$mapped_reads")

            # compute peak genome coverage
            peak_len=$(awk '{{total+=$3-$2}}END{{print total}}' $PEAK)
            genome_size={params.genomeSize}
            genomecov=$(bc -l <<< "$peak_len/$genome_size")

            # rounding fractions
            genomecov_round=$(printf "%.5f\n" "$genomecov")
            frip_round=$(printf "%.3f\n" "$frip")

            # write peak-based QC metrics to output file
            printf '{params.target}\\t'$CALLING_MODE'\\t'$peak_count'\\t'$frip_round'\\t'$genomecov_round'\\n' > {output.qc}

	        # add chr to peak files
            bedtools subtract -nonamecheck -a $PEAK -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > $PEAK_CHR
	        """

    rule aggregate_FRiP_SE:
        input:
            qc = expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAM_peaks.qc.txt", target = TARGETNAMES)
        output:
            aggregated_qc = "05_Quality_controls_and_statistics/peaks_stats/all_samples_FRiP_report.tsv"
        params:
            all_qc = ' '.join(expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAM_peaks.qc.txt", target = TARGETNAMES))
        threads: 1
        shell:
            """
            # print header of FRiP report
            printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
            cat {params.all_qc} >> {output.aggregated_qc}
            """


# ------------------------------------------------------------------------------
#                                 END pipeline
# ------------------------------------------------------------------------------
