function [drug_vector, apoptosis_ratio] = boolean_network_function(num_in,num_fault,num_drugs,inputs,faults,drugs)
%% Inputs : The input vector for the network is defined as [ER stress, TNF alpha, TRAIL, PTP, IL6, DNA Damage, IGF, EGF]

dec_in = inputs;

bin_in = dec2bin(dec_in,num_in);

for k = 1 : num_in
    input(k) = bin2dec(bin_in(k));
end

er_stress = input(1);
tnf_a = input(2);
trail = input(3);
ptp = input(4);
il6 = input(5);
dna_damage = input(6);
igf = input(7);
egf = input(8);

%% Faults : The fault vector is [Ras, Raf, PTEN, p53, STAT3, DR5]

for i = 1 : length(faults)
    
    dec_f = faults(i);
    bin_f = dec2bin(dec_f,num_fault);
    
    for k = 1 : num_fault
        fault(k) = bin2dec(bin_f(k));
    end
    
    faulty(i) = string(bin_f);
    
    ras_fault = fault(1);
    raf_fault = fault(2);
    pten_fault = fault(3);
    p53_fault = fault(4);
    stat3_fault = fault(5);
    dr5_fault = fault(6);
    
    stuck_at_0 = 0;
    stuck_at_1 = 1;
    
    %% Drugs : The drug vector consists of [Cryptotanshinone, LY294002, Temsirolimus, UO-126, Lapatinib, SH5-07, AG1024]
    
    for j = 1 : length(drugs)
        
        dec_dr = drugs(j);
        bin_dr = dec2bin(dec_dr,num_drugs);
        
        for k = 1 : num_drugs
            drug(k) = bin2dec(bin_dr(k));
        end
        
        drugy(j) = string(bin_dr);
        
        cryp = drug(1);
        ly294002 = drug(2);
        temsi =drug(3);
        uo126 = drug(4);
        lapatinib = drug(5);
        sh507 = drug(6);
        ag1024 = drug(7);
        
        %% PI3k-Akt Pathway
        
        igf1r = igf & (~ag1024);
        
        irs = igf1r;
        
        %% MAPK-ERK Pathway
        
        egfr = egf & (~lapatinib);
        
        if ras_fault == 0
            ras = egfr | irs;
        else
            ras = stuck_at_1;
        end
        
        if raf_fault == 0
            raf = ras;
        else
            raf = stuck_at_1;
        end
        
        mek1 = raf & (~uo126);
        erk1 = mek1;
        
        cmyc = erk1;
        cfos = erk1;
        
        mekk1 = ras;
        mkk7 = mekk1 | er_stress;
        ask1 = trail | tnf_a;
        jnk = mkk7 | ask1 ;
        
        pi3k = (ras | irs) & (~ly294002);
        pip2 = pi3k;
        
        if pten_fault == 0
            pten = 1; 
         else
            pten = stuck_at_0;
        end
        
        pip3 = pip2 & (~pten); 
        mtorc2 = pip3;
        pdk1 = pip3;
        akt = pdk1 | mtorc2;
        mtorc1 = akt & (~temsi);
        
        %% ER Stress Pathway
        
        casp12 = er_stress;
        
        perk = er_stress;
        
        atf4 = perk;
        
        elf2 = perk;
        
        %% STAT3 Pathway
        
        il6r = il6;
        jak1 = il6r | perk;
        
        temp = (~ptp)& (jak1 | erk1 | mtorc1);
        
        if stat3_fault == 0
            stat3 =   (temp & (~cryp))|(temp &(~sh507)); 
        else
            stat3 = (stuck_at_1 & (~cryp))|(stuck_at_1 &(~sh507)); 
        end
        
        mcl1 = stat3 | (~dna_damage);
        
        rad3 = dna_damage;
        
        vegf = stat3;
        
        chop = (atf4 | jnk)  & (~stat3);
        
        %% Cell Proliferation
        
        chk1 = mcl1 & rad3;
        
        mdm2 = erk1 | akt;
        
        if p53_fault == 0
            p53 = chk1 | (~mdm2) | (~jnk); 
        else
            p53 = stuck_at_0; 
        end
        
        %% p53 Pathway
        
        p21 = p53 & (~mtorc1);
        
        cyclind1_cdk4 = (~p21) & (mdm2 | cmyc | cfos);
        
        noxa = p53;
        
        puma = p53;
        
        %% NFKB Pathway
        
        rip = trail;
        
        tnfr = tnf_a;
        
        nfkb = tnfr | rip | elf2;
        
        cox2 = nfkb | stat3; 
        
        a1 = nfkb;
        
        %% TRAIL Apoptosis
        
        srp = trail;
        
        dr4 = srp & (~cox2);
        
        if dr5_fault == 0
            
            dr5 = (trail | (chop & (~cox2)) )| cryp ; 
        else
            dr5 = stuck_at_0 | cryp; 
        end
        
        disc = (dr4 | dr5) | tnfr;
        
        il8 = disc & nfkb & (~jnk); 
        
        cflip = il8;

       casp8 = mean(disc+(~cflip));
        
        %% Mitochondrial Apoptosis
        
        bid = casp8;
        
        bad = (~akt);
        
        bim = mean((~mcl1)+chop);
        
        bclxl = mean((~bad)+(~bid));
        
        bcl2 = mean((~bim)+(~noxa)+(~puma)+(puma & mcl1));
        
        bakx = mean((~bcl2)+(~bclxl)+bid+jnk);
        
        ros = mean(stat3+trail+tnf_a+er_stress);
        
        mmp = mean(bakx+ros+stat3+(~bcl2)+(~bclxl)+(~a1));
        
        cytoc = mmp;
        
        xiap = mean(akt + (~mmp));
        
        casp3 = (casp12 | cytoc) & (~xiap);
                
        %% Outputs
        
        angiogenesis(i,j) = mean(vegf + cox2);  
        
        cdc25 = (~chk1);
        
        wee1 = chk1;
        
        cdc2 = mean(cdc25 +(~wee1));
        
        cell_proliferation(i,j) = mean(cyclind1_cdk4+cdc2); 
        
        pro = mean(bakx + bad + casp8 + casp12 + bid + bim);
        anti = mean(bcl2 + bclxl + mcl1 + xiap);
        
        % The apoptosis ratio
        apop_ratio(i,j) = pro/anti;
        
    end
end

drug_vector = drugy;
apoptosis_ratio = apop_ratio;
end