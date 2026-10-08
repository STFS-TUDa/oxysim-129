/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::radiation::wsggGreyCloudAbsorptionEmission

Description
    Extension of the WSGG absorption emission model that accounts for
    radiation from Lagrangian carbonaceous particle clouds.

    The model:
      - Computes the volumetric particle emissive power from parcel
        temperatures, areas and a constant particle emissivity, and
        distributes it into the WSGG bands using gas weighting
        factors evaluated at the particle temperature and local
        mixture ratio.
      - Forms a gray volumetric particle absorption coefficient from
        the time-averaged particle surface area in each cell.

    The particle contributions are combined with the gas-phase
    radiation field via wsggBinaryAbsorptionEmission. Usage requires an
    enabled radiation submodel in the cloud properties dict, with an
    appropriate gas-phase WSGG coefficient set (e.g. Dorigon or Bordbar)
    in radiationProperties:

    absorptionEmissionModel wsggBinaryAbsorptionEmission;
    wsggBinaryAbsorptionEmissionCoeffs
    {
        model1
        {
            absorptionEmissionModel wsggDorigonAbsorptionEmission;

            wsggDorigonAbsorptionEmissionCoeffs
            {
                ...
            }
        }
        model2
        {
            absorptionEmissionModel wsggGreyCloudAbsorptionEmission;
            wsggGreyCloudAbsorptionEmissionCoeffs
            {
                cloudNames
                (
                    coalCloud
                );
            }
        }
    }

SourceFiles
    wsggGreyCloudAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_wsggGreyCloudAbsorptionEmission_H
#define radiation_wsggGreyCloudAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "fluidThermo.H"
#include "wsggAbsorptionEmission.H"
#include "radiation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggGreyCloudAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class wsggGreyCloudAbsorptionEmission
: public wsggAbsorptionEmission
{
    // Private data

    //- Coefficients dictionary
    dictionary coeffsDict_;

    //- Cloud name(s)
    const wordList cloudNames_;

    label nBands_;

    mutable PtrList<volScalarField> EDispfields;
    mutable PtrList<volScalarField> eDispfields;

    constexpr static double WCO2_ = 44.01;   // kg/kmol
    constexpr static double WH2O_ = 18.02;   // kg/kmol
    constexpr static double R_ = 8.314472e3; // ideal gas constant in m3*Pa/kmol/K

public:
    //- Runtime type information
    TypeName("wsggGreyCloudAbsorptionEmission");


    // Constructors

    //- Construct from components
    wsggGreyCloudAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggGreyCloudAbsorptionEmission();


    // Member Functions

    // Access

    // Absorption coefficient

    //- Absorption coefficient for dispersed phase
    virtual tmp<volScalarField> aDisp(const label bandI = 0) const;


    // Emission coefficient

    //- Emission coefficient for dispersed phase
    virtual tmp<volScalarField> eDisp(const label bandI = 0) const;


    // Emission contribution

    //- Return emission contribution for dispersed phase
    virtual tmp<volScalarField> EDisp(const label bandI = 0) const;

    // Emission coefficient

    //- Dispersed phase's contribution on continuos phase. It is 0. Required anyway due to unit overwritting
    virtual tmp<volScalarField> eCont(const label bandI = 0) const;

    // Number of bands --> Grey gases + 1 clear gas
    inline label nBands() const
    {
        return nBands_;
    }

    void correct(
        Foam::volScalarField& a,
        PtrList<Foam::volScalarField>& aLambda) const;

    template<class T>
    void particleSourceHelper(const T& tc) const;

    virtual Foam::scalarList getWeights(
        Foam::scalar T,
        Foam::scalar Yco2h2o) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
