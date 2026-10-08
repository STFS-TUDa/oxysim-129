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
    Foam::radiation::wsggDorigonAbsorptionEmission

Description
    Weighted-sum-of-gray-gases-model (WSGGM) for gas phase absorption/emission. Dorigon's model for air atmospheres.

Usage
    \verbatim
        wsggDorigonAbsorptionEmissionCoeffs //Coeffs and model reference Dorigon: http://dx.doi.org/10.1016/j.ijheatmasstransfer.2013.05.010
        {
            biCoeffs
            4
            (
                    5(0.5617e-1     7.844e-4        -8.563e-7       4.246e-10       -7.44e-14)
                    5(1.426e-1      1.795e-4        -0.1077e-7      -0.6972e-10     1.774e-14)
                    5(1.362e-1      2.574e-4        -3.711e-7       1.575e-10       -2.267e-14)
                    5(1.222e-1      -0.2327e-4      -0.74927e-7     0.4275e-10      -0.6608e-14)
            );

            kiCoeffs
            5
            (
                    0
                    1.921E-01
                    1.719E+00
                    1.137E+01
                    1.110E+02
            );
        }
    \endverbatim

SourceFiles
    wsggDorigonAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef wsggDorigonAbsorptionEmission_H
#define wsggDorigonAbsorptionEmission_H

#include "wsggAbsorptionEmission.H"
#include "fluidThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggDorigonAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class wsggDorigonAbsorptionEmission
: public wsggAbsorptionEmission
{

    // Private data

    //- Absorption model dictionary
    dictionary coeffsDict_;

    //- absorption coefficitent for individual grey gases
    scalarListList bi_;

    //- fitting Factors
    scalarList ki_;


    // characteristic path length
    // dimensionedScalar pathLength_;

    // Flag for meanBeamPath calculation
    // bool meanBeamPathAutoCalcMode_;

    // Domain sector //[degree]
    // dimensionedScalar sector_;
    constexpr static double WCO2_ = 44.01;   // kg/kmol
    constexpr static double WH2O_ = 18.02;   // kg/kmol
    constexpr static double R_ = 8.314472e3; // ideal gas constant in m3*Pa/kmol/K

    double Mr_;

    label nBands_;
    mutable PtrList<volScalarField> efields;
    mutable PtrList<volScalarField> afields;

public:
    //- Runtime type information
    TypeName("wsggDorigonAbsorptionEmission");


    // Constructors

    //- Construct from components
    wsggDorigonAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggDorigonAbsorptionEmission();


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
