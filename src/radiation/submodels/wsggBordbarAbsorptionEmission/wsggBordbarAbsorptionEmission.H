/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::radiation::wsggBordbarAbsorptionEmission

Description
    Weighted-sum-of-gray-gases-model (WSGGM) for gas phase absorption/emission. Bordbar's model for oxyfuel atmospheres.

Usage
    \verbatim
        wsggBordbarAbsorptionEmissionCoeffs //Coeffs and model reference Bordbar: http://dx.doi.org/10.1016/j.combustflame.2014.03.013
        {
            c1Coeffs
            5
            (
                    5(0.7412956     -0.5244441      0.582286        -0.2096994      0.0242031)
                    5(-0.9412652    0.2799577       -0.7672319      0.3204027       -0.0391017)
                    5(0.8531866     0.0823075       0.528943        -0.2468463      0.031094)
                    5(-0.3342806    0.1474987       -0.4160689      0.1697627       -0.0204066)
                    5(0.0431436     -0.0688622      0.1109773       -0.0420861      0.0049188)
            );

            c2Coeffs
            5
            (
                    5(0.1552073     -0.4862117      0.3668088       -0.1055508      0.0105857)
                    5(0.6755648     1.409271        -1.3834490      0.457521        -0.0501976)
                    5(-1.1253940    -0.5913199      0.9085441       -0.3334201      0.0384236)
                    5(0.6040543     -0.0553385      -0.1733014      0.0791608       -0.0098934)
                    5(-0.1105453    0.0464663       -0.0016129      -0.0035398      0.0006121)
            );

            c3Coeffs
            5
            (
                    5(0.2550242     0.3805403       -0.4249709      0.1429446       -0.0157408)
                    5(-0.6065428    0.3494024       0.1853509       -0.1013694      0.0130244)
                    5(0.8123855     -1.1020090      0.4046178       -0.0811822      0.0062981)
                    5(-0.4532290    0.6784475       -0.3432603      0.0883088       -0.0084152)
                    5(0.0869309     -0.1306996      0.0741446       -0.0202929      0.002011)
            );

            c4Coeffs
            5
            (
                    5(-0.0345199    0.2656726       -0.1225365      0.0300151       -0.0028205)
                    5(0.4112046     -0.5728350      0.292449        -0.0798076      0.0079966)
                    5(-0.5055995    0.4579559       -0.2616436      0.0764841       -0.0079084)
                    5(0.2317509     -0.1656759      0.1052608       -0.0321935      0.003387)
                    5(-0.0375491    0.022952        -0.0160047      0.0050463       -0.0005364)
            );

            dCoeffs
            4
            (
                    5(0.0340429     0.0652305       -0.0463685      0.0138684       -0.0014450)
                    5(0.3509457     0.7465138       -0.5293090      0.1594423       -0.0166326)
                    5(4.57074       2.168067        -1.4989010      0.4917165       -0.0542999)
                    5(109.8169      -50.923590      23.43236        -5.1638920      0.4393889)
            );
        }
    \endverbatim

SourceFiles
    wsggBordbarAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef wsggBordbarAbsorptionEmission_H
#define wsggBordbarAbsorptionEmission_H

#include "wsggAbsorptionEmission.H"
#include "fluidThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggBordbarAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class wsggBordbarAbsorptionEmission
: public wsggAbsorptionEmission
{

    // Private data

    //- Absorption model dictionary
    dictionary coeffsDict_;

    //- absorption coefficitent for individual grey gases
    scalarListList d_;

    //- fitting Factors
    scalarListList c1_;
    scalarListList c2_;
    scalarListList c3_;
    scalarListList c4_;

    const fvMesh& mesh_;

    // characteristic path length
    // dimensionedScalar pathLength_;

    // Flag for meanBeamPath calculation
    // bool meanBeamPathAutoCalcMode_;

    // Domain sector //[degree]
    // dimensionedScalar sector_;

    constexpr static double WCO2_ = 44.01;   // kg/kmol
    constexpr static double WH2O_ = 18.02;   // kg/kmol
    constexpr static double R_ = 8.314472e3; // ideal gas constant in m3*Pa/kmol/K


    DynamicList<scalarListList> Cijk;

    label nBands_;
    mutable PtrList<volScalarField> efields;
    mutable PtrList<volScalarField> afields;

public:
    //- Runtime type information
    TypeName("wsggBordbarAbsorptionEmission");


    // Constructors

    //- Construct from components
    wsggBordbarAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggBordbarAbsorptionEmission();


    // Member Functions

    // Access

    // Absorption coefficient

    //- Absorption coefficient for continuous phase
    tmp<volScalarField> aCont(const label bandI = 0) const;

    // Emission coefficient

    //- Emission coefficient for continuous phase
    tmp<volScalarField> eCont(const label bandI = 0) const;


    // Emission contribution

    //- Emission contribution for continuous phase
    tmp<volScalarField> ECont(const label bandI = 0) const;


    // Member Functions

    inline bool isGrey() const
    {
        return false;
    }

    // Number of bands --> Grey gases + 1 clear gas
    inline label nBands() const
    {
        return nBands_;
    }

    void correct(
        Foam::volScalarField& a,
        PtrList<Foam::volScalarField>& aLambda) const;

    Foam::scalarList getWeights(
        Foam::scalar T,
        Foam::scalar Xr) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
