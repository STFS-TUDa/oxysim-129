/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

inline const Foam::radiation::radiativeIntensityRaywsgg&
Foam::radiation::fvDOMwsgg::IRay(const label rayI) const
{
    return IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiation::fvDOMwsgg::IRayLambda(
    const label rayI,
    const label lambdaI) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiation::fvDOMwsgg::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiation::fvDOMwsgg::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiation::fvDOMwsgg::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiation::fvDOMwsgg::nLambda() const
{
    return nLambda_;
}


inline Foam::label Foam::radiation::fvDOMwsgg::nBands() const
{
    return nLambda_;
}

inline const Foam::volScalarField& Foam::radiation::fvDOMwsgg::a() const
{
    return a_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOMwsgg::aLambda(
    const label lambdaI) const
{
    return aLambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiation::fvDOMwsgg::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOMwsgg::qr() const
{
    return qr_;
}

inline const Foam::volScalarField& Foam::radiation::fvDOMwsgg::qin() const
{
    return qin_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOMwsgg::qem() const
{
    return qem_;
}


inline const Foam::radiation::blackBodyEmission&
Foam::radiation::fvDOMwsgg::blackBody() const
{
    return blackBody_;
}


inline Foam::scalar Foam::radiation::fvDOMwsgg::omegaMax() const
{
    return omegaMax_;
}


inline Foam::vector Foam::radiation::fvDOMwsgg::meshOrientation() const
{
    return meshOrientation_;
}


inline bool Foam::radiation::fvDOMwsgg::useSolarLoad() const
{
    // Solar load support (solarLoad) was stripped from this port — see the
    // note in fvDOMwsgg.H. Always false so callers take the no-solar-load path.
    return false;
}


inline bool Foam::radiation::fvDOMwsgg::useExternalBeam() const
{
    return useExternalBeam_;
}


inline const Foam::scalarList& Foam::radiation::fvDOMwsgg::
    spectralDistribution() const
{
    return spectralDistribution_;
}


// ************************************************************************* //
