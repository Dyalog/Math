:Namespace complex

    ⎕ML←1                         ⍝ ∊ is enlist
    Join←{⍺⍪⍉n n⍴⍵}               ⍝ hang vectors below values
    Real←1289≠⎕dr                 ⍝ is mat real?
    Sym←{w≡+⍉w←⍵+(○+*)≢⎕CT←2*¯32} ⍝ offset coarse tolerance for max fuzz symmetric/hermetian check.
    J←{⍺+¯11○⍵}                   ⍝ J's j. verb
    V2J←J/¨                       ⍝ convert old vector notation (matrix of real-cmpx pairs) to new 1J2 notation
    J2V←9 11∘○¨                   ⍝ convert new 1J2 notation to old vector notation (matrix of real-cmpx pairs)
    Call←{(⍎⍺⍺ Assoc ⍵⍵)⍵}        ⍝ Call external fn ⍺⍺ with arg ⍺, associate with types ⍵⍵ if necessary
    Assoc←{3=⎕NC ⍺:⍺ ⋄ ⎕NA'lapack',(¯2↑'32',⎕D∩⍨⊃# ⎕WG'APLVersion'),'|',⍺,' ',⍵} ⍝ associate external function.

      Dsyev←{
    ⍝  JOBZ  UPLO  N  A  LDA  W  WORK  LWORK  INFO
          types←'<C1 <C1 <I4 =F8[] <I4 >F8[] >F8[] <I4 >I4'
          args←'V' 'L'n(,⍵)n n(¯1+3×n)(¯1+3×n)0
          2 1 3⊃¨⊂1 1 0 1/('dsyev_'Call types)args ⍝ call external fn.
      }

      Dgeev←{
   ⍝  JOBVL  JOBVR  N  A  LDA  WR  WI  VL  LDVL  VR  LDVR  WORK  LWORK  INFO
          types←'<C1 <C1 <I4 =F8[] <I4 >F8[] >F8[] >I4 <I4 >F8[] <I4 >F8[] <I4 >I4'
          args←'N' 'V'n(,⍉mat)n n n 0 1(n×n)n(4×n)(4×n)0
          (real cmpx vecs code)←0 1 1 0 1 0 1/('dgeev_'Call types)args ⍝ call external fn.
     
          vals←real J cmpx        ⍝ combine parts
          ∧/0=cmpx:real vecs code ⍝ all-real case
          vr←⍉n n⍴vecs            ⍝ build matrix
          eii←vr×⍤1⊢cmpx=0
          wi←cmpx,0                           ⍝ {(⍵≠0)^⍵=-1⌽⍵}WI,0
          pair←(wi≠0)∧wi=-1⌽wi                ⍝ problem if 3 same in a row!
          pair←pair/⍳⍴pair
          eii[;pair]←vr[;pair]J vr[;1+pair]
          eii[;1+pair]←vr[;pair]J-vr[;1+pair]
          vecs←⍉eii
          vals vecs code                      ⍝ return parts
      }

      Zheev←{
       ⍝  JOBZ  UPLO  N  A  LDA  W  WORK  LWORK  RWORK  INFO
          types←'<C1 <C1 <I4 =J16[] <I4 >F8[] >F8[] <I4 >F8[] >I4'
          args←'V' 'L'n(∊⍉mat)n n(¯2+4×n)(¯1+2×n)(¯2+3×n)0
          2 1 3⊃¨⊂1 1 0 0 1/('zheev_'Call types)args ⍝ call external fn.
      }

      Zgeev←{
       ⍝  JOBVL  JOBVR  N  A  LDA  W  VL  LDVL  VR  LDVR  WORK  LWORK  RWORK  INFO
          types←'<C1 <C1 <I4 =J16[] <I4 >J16[] <I4 <I4 >J16[] <I4 >F8[] <I4 >F8[] >I4'
          args←'N' 'V'n(∊⍉mat)n n 0 1(n×n)n(4×n)(2×n)(2×n)0
          0 1 1 0 0 1/('zgeev_'Call types)args ⍝ call external fn.
      }

    ∇ vv←Eigen mat;vecs;vals;code;n
      n←≢mat
     
      :If Real mat
          :If (Sym mat)∧2≠n ⍝ symmetric
              (vals vecs code)←Dsyev mat
     
          :Else ⍝ general
              (vals vecs code)←Dgeev mat
          :EndIf
     
      :Else ⍝ complex
          :If Sym mat ⍝ hermitian
              (vals vecs code)←Zheev mat
     
          :Else ⍝ general
              (vals vecs code)←Zgeev mat
          :EndIf
     
      :EndIf
      ⎕SIGNAL code/11
      vv←vals Join vecs
    ∇

    Fourier←{⍺←1 ⋄ ⍵ RunFn 'idft' 'dft'⊃⍨1+0⌈⍺}
    ∇ data←data RunFn func;rank;shape;count ⍝ n-D Discrete Fourier Transformation
      (rank count)←(⍴,×/)shape←⍴data
      :If 1<count ⍝ non-scalar
          :If 3≠⎕NC func ⍝ associate external function if it does not exist
              ⎕NA'fftw',(¯2↑'32',⎕D∩⍨⊃# ⎕WG'APLVersion'),'|',func,'<I4 <I4[] =J16[]'
          :EndIf
          data←(count*0.5)÷⍨shape⍴(⍎func)rank shape,⊂∊data ⍝ Normalized Fourier/Inverse Fourier transform
      :EndIf
    ∇

    :namespace test

        ∇ ok←QA
          ok←eigen∧fourier∧bfour
        ∇

        symm←{⍵+⍉⍵}             ⍝ symmetrify
        herm←⍉++                ⍝ hermitianise
        rand←{(-⍵×5)+?⍵ ⍵⍴⍵×10} ⍝ random matrix size ⍵
        posdef←{s←≢⍵ ⋄ ↑+.×/(⍉⍵)(s s⍴↑(s+1)↑¨0.1×s?100)⍵} ⍝ positive definite
        comp←{⍵+¯11○rand≢⍵}     ⍝ complexify

        ∇ ok←bfour;mats;mat;Test;complex
        ⍝ runs 50 matrices checking that everything is as before
          ok←1
          Test←##.Fourier≡##.V2J∘#.Fourier∘##.J2V
          mats←rand¨⍳50
          :For mat :In mats
              ok∧←Test symm mat             ⍝ real symmetric
              ok∧←Test posdef mat           ⍝ real positive definite
              ok∧←Test mat                  ⍝ real
              ok∧←Test herm comp mat        ⍝ complex hermitian
              ok∧←Test posdef herm comp mat ⍝ complex positive definite
              ok∧←Test comp mat             ⍝ complex
          :EndFor
        ∇

        ∇ ok←fourier;⎕PATH;⎕IO;x;fuzz;th;n;Xn;Yn;Zn;Kn;Ln;In;Inn;Enn;i;j  ⍝ Complex 1-D test of Fourier
        ⍝ runs 50 random matrices
          ok←1
          ⎕PATH,←' ↑'
          ⎕IO←1
          x←{                         ⍝ complex multiply
              sign←2 2⍴1 ¯1 1 1
              +/sign×0 1⊖(2↑⍺)∘.×2↑⍵
          }
          fuzz←{⎕CT←2*¯32 ⋄ (⍺+○1)⍺⍺ ⍵+○1}
          :For i :In 10?200
              :For j :In 5?100
                  th←(⍳i)÷j ⍝ th←(⍳?200)÷?100
                  ⍝'' ⋄ 'n is ',⍕⍴th
                  Xn←2↑¨⊃+/,1 2∘.○(⊂th)×¨6?6 ⍝ try a step function like Xn←2↑¨(50⍴1),(50⍴2)
                  ⍝------------------------------
                  Kn←##.J2V ##.Fourier ##.V2J Xn
                  ⍝now test answer
                  n←⍴Xn
                  In←(⍳n)-⎕IO
                  ⍝Enn←*○0J¯2×In∘.×In÷n
                  Inn←○2×In∘.×In÷n
                  Enn←↓(2○Inn),[⎕IO+1.5](-1○Inn)
                  Ln←(÷0.5*⍨n)×Xn+.x Enn
                  ⍝'Does APL model give same answers?'
                  ok∧←∧/∊Ln=fuzz Kn ⍝ 'APL≡Fourier?'
                  ⍝------------------------------
                  Yn←##.J2V ¯1 ##.Fourier ##.V2J Kn
                  ⍝now test answer
                  ⍝'Does Inverse Fourier return original argument?'
                  ok∧←∧/∊Xn=fuzz Yn
                  ⍝Enn←*○0J2×In∘.×In÷n
                  Enn←↓(2○Inn),[⎕IO+1.5](1○Inn)
                  Zn←(÷0.5*⍨n)×Kn+.x Enn
                  ⍝'Does APL model of inverse give same answers?'
                  ok∧←∧/∊Zn=fuzz Yn ⍝ 'APL≡Inverse Fourier?'
                  ⍝------------------------------
              :EndFor
          :EndFor
        ∇

        ∇ ok←eigen;⎕CT;⎕PP;i;Aii;Eji;real;fuzz;⎕PATH;x;norm;conj;anti;symm;herm;diag;j;k;l
        ⍝ 480 random matrices
          ok←1
          ⎕PATH,←' ↑'
          x←{                    ⍝ complex multiply
              sign←2 2⍴1 ¯1 1 1
              +/sign×0 1⊖(2↑⍺)∘.×2↑⍵
          }
          real←{1=≡,⍵}           ⍝ matrix is real.
          conj←{⍵×⊂1 ¯1}         ⍝ complex conjugate.
          anti←{⍵-⍉⍵}            ⍝ antisymmetrify
          symm←{⍵+⍉⍵}            ⍝ symmetrify
          herm←{real cmpx←⊂[1 2]↑⍵ ⋄ (symm real),¨anti cmpx}  ⍝ hermitianise.
          :For i :In 10/⍳6
              :For j :In 0 1
                  :For k :In 0 1
                      :For l :In 0 1
                          :If j  ⍝ real
                              Aii←(-i÷2)+?i i⍴i               ⍝ real non-symmetric
                              :If k   ⍝ symmetric
                                  :If l        ⍝ not nec pos def
                                      ⍝ 'real symmetric'
                                      Aii+←⍉Aii                   ⍝ real symmetric
                                  :Else
                                      ⍝ 'real positive definite'
                                      diag←i i⍴↑(i+1)↑¨0.1×i?100
                                      Aii←↑+.×/(⍉Aii)diag Aii
                                  :End
                              :Else
                                  ⍝ 'real non-symmetric'
                              :End
                          :Else      ⍝ complex
                              Aii←↓(-i÷2)+?i i 2⍴i
                              :If k   ⍝ hermitian
                                  :If l    ⍝ not nec positive definite.
                                      ⍝ 'complex hermitian'
                                      Aii←herm Aii            ⍝ complex  hermitian
                                  :Else
                                      ⍝ 'complex positive definite'
                                      Aii←herm Aii
                                      diag←2↑¨i i⍴↑(i+1)↑¨0.1×i?100
                                      Aii←↑+.x/(⍉conj Aii)diag Aii     ⍝ complex positive definite
                                  :End
                              :Else
                                  ⍝ 'complex non-hermitian'
                              :End
                          :End
                          :Trap 11
⍝-------------------------------------------------------------
                              Eji←##.J2V ##.Eigen ##.V2J Aii       ⍝N.B. should choose?  Eii←Eii××+/×,Eii
⍝-------------------------------------------------------------
                          ⍝⌊0.5+⊃↓Eji ⍝ check pos def case
                          :Else
                              ok←0 ⋄ 'failed' ⋄ Aii ⋄ →0        ⍝try DSYEV→DGEEV/ZHEEV→ZGEEV
                          :End
                          norm←{                         ⍝ Euclidean norm
                              1=≡,⍵:(+⌿⍵*2)*0.5          ⍝ real
                              (+⌿,[1 2]2 3 1⍉↑⍵*2)*0.5   ⍝ complex
                          }
                          fuzz←{⎕CT←2*¯32 ⋄ (⍺+○1)⍺⍺ ⍵+○1}
         
                          :If real Eji
                              ⍝'norm'(norm 1 0↓Eji)
                              ok∧←∧/∊(Aii+.×1 0↓Eji)=fuzz(Eji[1;]×[2]1 0↓Eji)
                          :Else
                              :If real Aii
                                  Aii←2↑¨Aii
                              :End
                              ⍝'euclidean norm'(norm 1 0↓Eji)
                              ok∧←∧/∊(Aii+.x 1 0↓Eji)=fuzz(i i⍴Eji[1;])x¨1 0↓Eji
                          :End
                      :EndFor
                  :EndFor
              :EndFor
          :EndFor
        ∇
    :endnamespace

:EndNamespace