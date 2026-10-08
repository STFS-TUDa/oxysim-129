/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "C2SM.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::C2SM<CloudType>::
    C2SM(
        const dictionary& dict,
        CloudType& owner)
: DevolatilisationModelSTFS<CloudType>(dict, owner, typeName),
  volatileData_(this->coeffDict().lookup("volatileData")),
  YVolatile0_(volatileData_.size()),
  volatileToGasMap_(volatileData_.size()),
  residualCoeff_(this->coeffDict().getScalar("residualCoeff"))
{
    if (volatileData_.empty())
    {
        WarningInFunction
            << "Devolatilisation model selected, but no volatiles defined"
            << nl << endl;
    }
    else
    {
        Info << "Participating volatile species:" << endl;

        // Determine mapping between active volatiles and cloud gas components
        const label idGas = owner.composition().idGas();
        const scalar YGasTot = owner.composition().YMixture0()[idGas];
        const scalarField& YGas = owner.composition().Y0(idGas);
        forAll(volatileData_, i)
        {
            const word& specieName = volatileData_[i].name();
            const label id = owner.composition().localId(idGas, specieName);
            volatileToGasMap_[i] = id;
            YVolatile0_[i] = YGasTot * YGas[id];

            Info << "    " << specieName << ": particle mass fraction = "
                 << YVolatile0_[i] << endl;

            Info << "Gas constant defined as: " << RR << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::C2SM<CloudType>::calculate(
    const scalar dt,
    const scalar age,
    const scalar mass0,
    const scalar mass,
    const scalar T,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV,
    scalarField& dMassC) const
{
    // bool done = true;
    forAll(volatileData_, i)
    {
        const label id = volatileToGasMap_[i];
        const scalar massVolatile0 = mass0 * YVolatile0_[i];
        const scalar massVolatile = mass * YGasEff[id];

        // Combustion allowed once all volatile components evolved
        // done = done && (massVolatile <= residualCoeff_*massVolatile0);

        // Model coefficients

        const scalar a1 = volatileData_[i].a1();
        const scalar A11 = volatileData_[i].A11();
        const scalar E1 = volatileData_[i].E1();
        const scalar a2 = volatileData_[i].a2();
        const scalar A12 = volatileData_[i].A12();
        const scalar E2 = volatileData_[i].E2();

        const scalar kappa1 = A11 * exp(-E1 / (RR * T));
        const scalar kappa2 = A12 * exp(-E2 / (RR * T));

        // Mass transferred from particle to carrier gas phase
        dMassDV[id] = min(dt * massVolatile * (kappa1 * a1 + kappa2 * a2), (kappa1 * a1 + kappa2 * a2) / (kappa1 + kappa2) * massVolatile);
        dMassC[id] = min(dt * massVolatile * (kappa1 * (1 - a1) + kappa2 * (1 - a2)), massVolatile - dMassDV[id]);
    }

    // if (done && canCombust != -1)
    //{
    canCombust = 1;
    //}
}


// ************************************************************************* //
