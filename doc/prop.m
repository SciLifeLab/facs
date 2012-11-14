hit =[1670777
464386
131302
44158
22436
17064
]

hit1 = [10342350
3092821
842175
241199
88697
50011
]

unhit = [47329223
48335614
48468698
48355842
48177564
47982936
]

unhit1 = [39457650
46507179
48557825
48958801
48911303
48749989
]

total = hit+unhit
total1 = hit1+unhit1

polp = hit./total
polp1 = hit1./total1
hold on
plot(polp,'-rs');
plot(polp1,'-ks');
hold off