function integrandFunct = getIntegrand(iFunct1, iFunct2, iFunct3, jFunct1, jFunct2, jFunct3, basisICoeffs, basisJCoeffs, r, z, k, idNum)


    curlRZKforPhiLine1 = @(bCoeff) [   bCoeff   ];
    curlRZKforPhiLine2 = @(bCoeff) [ - bCoeff/k ];
    curlRZKforPhiLine3 = @(aCoeff) [   aCoeff   ];
                      
    curlRZKforPsiLine1 = @(coeffBasis) [ - coeffBasis(3) - coeffBasis(1).*r          ];
    curlRZKforPsiLine2 = @(coeffBasis) [ - coeffBasis(3)/k - 3.*(coeffBasis(1)/k).*r ];
    curlRZKforPsiLine3 = @(coeffBasis) [   coeffBasis(2) - coeffBasis(1).*z          ];
    
    
    
    curlsDotProductFor1 = @() [  ( curlRZKforPhiLine1(basisICoeffs(2)).*curlRZKforPhiLine1(basisJCoeffs(2)) ) + ( curlRZKforPhiLine2(basisICoeffs(2)).*curlRZKforPhiLine2(basisJCoeffs(2)) ) + ( curlRZKforPhiLine3(basisICoeffs(1)).*curlRZKforPhiLine3(basisJCoeffs(1)) )  ];
    curlsDotProductFor2 = @() [  ( curlRZKforPhiLine1(basisICoeffs(2)).*curlRZKforPsiLine1(basisJCoeffs)    ) + ( curlRZKforPhiLine2(basisICoeffs(2)).*curlRZKforPsiLine2(basisJCoeffs)    ) + ( curlRZKforPhiLine3(basisICoeffs(1)).*curlRZKforPsiLine3(basisJCoeffs)    )  ];
    curlsDotProductFor3 = @() [  ( curlRZKforPsiLine1(basisICoeffs).*curlRZKforPhiLine1(basisJCoeffs(2))    ) + ( curlRZKforPsiLine2(basisICoeffs).*curlRZKforPhiLine2(basisJCoeffs(2))    ) + ( curlRZKforPsiLine3(basisICoeffs).*curlRZKforPhiLine3(basisJCoeffs(1))    )  ];
    curlsDotProductFor4 = @() [  ( curlRZKforPsiLine1(basisICoeffs).*curlRZKforPsiLine1(basisJCoeffs)       ) + ( curlRZKforPsiLine2(basisICoeffs).*curlRZKforPsiLine2(basisJCoeffs)       ) + ( curlRZKforPsiLine3(basisICoeffs).*curlRZKforPsiLine3(basisJCoeffs)       )  ];
                             
                             
    
    
    if idNum == 1
        integrandFunct = [ ( curlsDotProductFor1() + (iFunct1.*jFunct1 + iFunct2.*jFunct2 + iFunct3.*jFunct3) ) .*r];
    end
    
    if idNum == 2
        integrandFunct = [ ( curlsDotProductFor2() + (iFunct1.*jFunct1 + iFunct2.*jFunct2 + iFunct3.*jFunct3) ) .*r];
    end
    
    if idNum == 3
        integrandFunct = [ ( curlsDotProductFor3() + (iFunct1.*jFunct1 + iFunct2.*jFunct2 + iFunct3.*jFunct3) ) .*r];
    end
    
    if idNum == 4
        integrandFunct = [ ( curlsDotProductFor4() + (iFunct1.*jFunct1 + iFunct2.*jFunct2 + iFunct3.*jFunct3) ) .*r];
    end
    
end