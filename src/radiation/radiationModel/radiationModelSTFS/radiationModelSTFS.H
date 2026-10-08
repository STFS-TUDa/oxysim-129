/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

Namespace
    Foam::radiation

Description
    Namespace for radiation modelling

Class
    Foam::radiation::radiationModelSTFS

Description
    Top level model for radiation modelling

SourceFiles
    radiationModelSTFS.C
    radiationModelSTFSNew.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelSTFS_H
#define radiationModelSTFS_H

#include "radiationModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatricesFwd.H"
#include "Switch.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;
class fvMesh;

namespace radiation
{

// Forward declaration of classes
// class absorptionEmissionModel;
class scatterModel;
class sootModel;

/*---------------------------------------------------------------------------*\
                       Class radiationModelSTFS Declaration
\*---------------------------------------------------------------------------*/

class radiationModelSTFS
: public radiationModel
{
public:
protected:
    //	volScalarField ShGas;


private:
    mutable DimensionedField<scalar, volMesh> ShGas;

    //- Pointer to the thermo model
    mutable const basicThermo* thermoPtr_;

    //- Resolve which basicThermo object radiation should use: an explicit
    //  'thermoName' entry in the given dict, else the sole registered
    //  basicThermo, else a fatal error naming the ambiguity/absence.
    static Foam::word resolveThermoName(
        const Foam::dictionary& dict,
        const Foam::fvMesh& mesh);

public:
    // Constructors

    //- Null constructor
    radiationModelSTFS(const volScalarField& T);

    //- Construct from components
    radiationModelSTFS(const word& type, const volScalarField& T);

    //- Construct from components
    radiationModelSTFS(
        const word& type,
        const dictionary& dict,
        const volScalarField& T);


    //- Destructor
    virtual ~radiationModelSTFS();


    // Member Functions

    //- Correct radiation model
    virtual void correct();

    // Edit

    //- Solve radiation equation(s)
    virtual void calculate() = 0;

    // Read radiationProperties dictionary is inherited as-is from
    // radiationModel (no longer redeclared pure-virtual here).

    //- Source term component (for power of T^4)
    virtual tmp<volScalarField> Rp() const = 0;

    //- Source term component (constant)
    virtual tmp<volScalarField::Internal> Ru() const = 0;

    //- Energy source term
    virtual tmp<fvScalarMatrix> Sh(
        const basicThermo& thermo,
        const volScalarField& he) const;

    virtual label nBands() const = 0;

    //- Access to the resolved thermo model. Fatal errors if none has been
    //  set (this should never happen post-construction).
    const basicThermo* thermoPtr() const
    {
        if (!thermoPtr_)
        {
            FatalErrorInFunction
                << "No thermo has been set on this radiationModelSTFS."
                << exit(FatalError);
        }
        return thermoPtr_;
    }
};

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
