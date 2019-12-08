"""Expected input: a file called 'exons.txt' containing four lines of sense DNA sequences.
Each line should correspond to the sequences for the exons 2,3,4, and 5 in that order.
The program will return the expected spectral sensitivity of the opsin encoded by the DNA
as well as a predicted classification (L or M)"""

tyrosine = ['TAT','TAC']
threonine = ['ACT','ACC', 'ACA', 'ACG']
phenylalanine = ['TTT', 'TTC']
alanine = [‘GCT’, ‘GCC’, ‘GCA’, ‘GCG’]
serine = [‘TCT’, ‘TCC’, ‘TCA’, ‘TCG’]
isoleucine = [‘ATT’, ’ATC’, ‘ATA’]
tyrosine = [‘TAT’, ‘TAC’]
opsinType = 'undefined'
baseSensitivity = 559 #nm sensitivity of wild-type L opsin

"""Helper function to return the codon at a specified position. 
Expected to be called in context of a given exon (string).
Sanity-check: the first codon of the exon should be the 'start' codon, 'ATG' for sense-DNA strand."""
def codonAt(exon, position):
    assert exon[0:3] == 'ATG', 'This is not a properly formatted sense-DNA exon.'
    return exon[3*(position-1):3*position]

"""Rules for exon 5: these determine overall L vs M type
L pigments have the amino acid tyrosine (‘TAT’ or ‘TAC’) at codon positions 277 and 309, 
and threonine (‘ACT’, ‘ACC’, ‘ACA’, or ‘ACG’) at codon position 285.
M pigments have phenylalanine (‘TTT’ or ‘TTC’) at codon positions 277 and 309, 
and alanine (‘GCT’, ‘GCC’, ‘GCA’, or ‘GCG’) at 285.
Substituting both the amino acid at 277 from tyrosine to phenylalanine 
and threonine to alanine at 285 results in a 20nm shortening of peak sensitivity; 
the literature is not clear on the individual contributions at either codon position.
Phenylalanine instead of tyrosine at codon 309 will shorten the peak sensitivity only 1-2nm."""
def ex5Rules(exon5,baseSensitivity):
    assert exon5[0:3] == 'ATG', 'This is not a properly formatted sense-DNA exon.'
    c277 = codonAt(exon5,277)
    c309 = codonAt(exon5, 309)
    c285 = codonAt(exon5, 285)
    if c277 in phenylalanine:
        if c285 in alanine:
            baseSensitivity-=20
            opsinType = 'M'
            print('Exon5 suggests M-type opsin')
        elif c285 in threonine:
            print('Indeterminate hybrid L/M type opsin; sensitivity unclear')
        else:
            print('Unexpected amino acid at codon position 285 of exon5'))
    elif c277 in tyrosine:
        if c285 in alanine:
            print('Indeterminate hybrid L/M type opsin; sensitivity unclear')
        elif c285 in threonine:
            #no change in base sensitivity
            opsinType = 'L'
            print('Exon5 suggests L-type opsin')
    else: 
        print('Unexpected amino acid at codon position 277 of exon5'))

    if c309 in phenylalanine:
        baseSensitivity-=2
        print('Exon 5 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
    elif c309 in tyrosine:
        #baseSensitivity unchanged
        print('Exon 5 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
    else:
        opsinType = 'undefined'
        print('Unexpected amino acid at codon position 309 of exon5')
    return 

"""Rules for exon 2: these further tune the sensitivity
Exon 2 of L opsins exhibits three polymorphisms at codon positions 65, 111, and 116, 
but only codon position 116 has an effect on spectral sensitivity.
L opsin DNA sequences with serine (‘TCT’, ‘TCC’, ‘TCA’, or ‘TCG’) 
at position at 116 result in sensitivities 2-3nm longer (red-shifted) 
than L pigments with tyrosine (‘TAT’ or ‘TAC’) at the same codon position. 
This substitution has negligible effects in M opsin sequences."""
def ex2Rules(exon2, opsinType, baseSensitivity):
    assert exon2[0:3] == 'ATG', 'This is not a properly formatted sense-DNA exon.'
    assert opsinType!='undefined', 'Rules for exon 2 only significantly apply to L-type opsins'
    if opsinType == 'M':
        #baseSensitivity unchanged
        print('Exon 2 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
    elif opsinType == 'L':
        c116 = codonAt(exon2,116)
        if c116 in serine:
            #baseSensitivity unchanged
            print('Exon 2 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
        elif c116 in tyrosine:
            baseSensitivity -=2
            print('Exon 2 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
        else:
            opsinType = 'undefined'
            print('Unexpected amino acid at codon position 116 of exon2')
    else:
        print('Something has gone horribly awry!')
    return

"""Rules for exon 3: these further tune the sensitivity
Exon 3 has the highest variability, but only the presence of serine or alanine at codon 180 
has an impact on spectral sensitiv- ity - serine results in a longer peak sensitivity wavelength. 
The substitution, similar to that of the serine/tyrosine substitution in exon 2, 
has a larger effect in L opsins than in M opsins."""
def ex3Rules(exon3, opsinType, baseSensitivity):
    assert exon3[0:3] == 'ATG', 'This is not a properly formatted sense-DNA exon.'
    assert opsinType!='undefined', 'Rules for exon 3 only significantly apply to L-type opsins'
    if opsinType == 'M':
        #baseSensitivity unchanged
        print('Exon 3 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
    elif opsinType == 'L':
        c180 = codonAt(exon3,180)
        if c180 in serine:
            #baseSensitivity unchanged
            print('Exon 3 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
        elif c180 in alanine:
            baseSensitivity -=2
            print('Exon 3 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
        else:
            opsinType = 'undefined'
            print('Unexpected amino acid at codon position 180 of exon3')
    else:
        print('Something has gone horribly awry!')
    return

"""Rules for exon 4: further tuning
Isoleucine (‘ATT’, ’ATC’, or ‘ATA’) at position 230 and alanine at position 233 
are red shifted compared to pigments with threonine at 230 and isoleucine at 233. 
As with exons 2 and 3, the shifts are larger for L opsins than M opsins, 
and the literature is not overly precise about the exact shifts."""
def ex4Rules(exon4, opsinType, baseSensitivity):
    assert exon4[0:3] == 'ATG', 'This is not a properly formatted sense-DNA exon.'
    assert opsinType!='undefined', 'Rules for exon 4 only significantly apply to L-type opsins'
    if opsinType == 'M':
        #baseSensitivity unchanged
        print('Exon 4 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
    elif opsinType == 'L':
        c230 = codonAt(exon4,230)
        c233 = codonAt(exon4,233)
        if c230 in isoleucine:
            if c233 in alanine:
                #baseSensitivity unchanged
                print('Exon 4 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
            elif c233 in isoleucine:
                print('Exon 4 suggests intermediate configuration, exact sensitivity unknown but roughly ' + baseSensitivity + 'nm')
            else:
                print('Unexpected amino acid at position 233 of exon 4.')
        elif c230 in threonine:
            if c233 in alanine:
                print('Exon 4 suggests intermediate configuration, exact sensitivity unknown but roughly ' + baseSensitivity + 'nm')
            elif c233 in isoleucine:
                baseSensitivity -=2
                print('Exon 4 suggests estimated peak spectral sensitivity is ' + baseSensitivity + 'nm')
            else:
                print('Unexpected amino acid at position 233 of exon 4.')
        else:
            print('Unexpected amino acid at position 230 of exon4')
    else:
        print('Something has gone horribly awry!')
    return


with open('exons.txt') as mysteryOpsin:
    exons = mysteryOpsin.readlines() #array object containing each of the lines
ex5Rules(exons[4],baseSensitivity)
ex2Rules(exons[0],opsinType,baseSensitivity)
ex3Rules(exons[1],opsinType,baseSensitivity)
ex4Rules(exons[2],opsinsType,baseSensitivity)