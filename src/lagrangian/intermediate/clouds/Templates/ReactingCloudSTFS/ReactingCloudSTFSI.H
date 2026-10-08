/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::ReactingCloudSTFS<CloudType>&
Foam::ReactingCloudSTFS<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::ReactingCloudSTFS<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::ReactingCloudSTFS<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::CompositionModel<Foam::ReactingCloudSTFS<CloudType>>&
Foam::ReactingCloudSTFS<CloudType>::composition() const
{
    return *compositionModel_;
}


template<class CloudType>
inline const Foam::PhaseChangeModel<Foam::ReactingCloudSTFS<CloudType>>&
Foam::ReactingCloudSTFS<CloudType>::phaseChange() const
{
    return *phaseChangeModel_;
}


template<class CloudType>
inline Foam::PhaseChangeModel<Foam::ReactingCloudSTFS<CloudType>>&
Foam::ReactingCloudSTFS<CloudType>::phaseChange()
{
    return *phaseChangeModel_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloudSTFS<CloudType>::rhoTrans(const label i)
{
    return rhoTrans_[i];
}


template<class CloudType>
inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::ReactingCloudSTFS<CloudType>::rhoTrans() const
{
    return rhoTrans_;
}


template<class CloudType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::ReactingCloudSTFS<CloudType>::rhoTrans()
{
    return rhoTrans_;
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix> Foam::ReactingCloudSTFS<CloudType>::SYi(
    const label i,
    volScalarField& Yi) const
{
    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("Yi"))
        {
            tmp<volScalarField> trhoTrans(
                new volScalarField(
                    IOobject(
                        this->name() + ":rhoTrans",
                        this->db().time().timeName(),
                        this->db(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false),
                    this->mesh(),
                    dimensionedScalar(dimMass / dimTime / dimVolume, Zero)));

            volScalarField& sourceField = trhoTrans.ref();

            sourceField.primitiveFieldRef() =
                rhoTrans_[i] / (this->db().time().deltaTValue() * this->mesh().V());

            const dimensionedScalar YiSMALL("YiSMALL", dimless, SMALL);

            return fvm::Sp(neg(sourceField) * sourceField / (Yi + YiSMALL), Yi)
                 + pos0(sourceField) * sourceField;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Yi, dimMass / dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            fvm.source() = -rhoTrans_[i] / this->db().time().deltaTValue();

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>::New(Yi, dimMass / dimTime);
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ReactingCloudSTFS<CloudType>::Srho(const label i) const
{
    tmp<volScalarField::Internal> tRhoi(
        new volScalarField::Internal(
            IOobject(
                this->name() + ":rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false),
            this->mesh(),
            dimensionedScalar(
                rhoTrans_[0].dimensions() / dimTime / dimVolume, Zero)));

    if (this->solution().coupled())
    {
        scalarField& rhoi = tRhoi.ref();
        rhoi = rhoTrans_[i] / (this->db().time().deltaTValue() * this->mesh().V());
    }

    return tRhoi;
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ReactingCloudSTFS<CloudType>::Srho() const
{
    tmp<volScalarField::Internal> trhoTrans(
        new volScalarField::Internal(
            IOobject(
                this->name() + ":rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false),
            this->mesh(),
            dimensionedScalar(
                rhoTrans_[0].dimensions() / dimTime / dimVolume, Zero)));

    if (this->solution().coupled())
    {
        scalarField& sourceField = trhoTrans.ref();
        forAll(rhoTrans_, i)
        {
            sourceField += rhoTrans_[i];
        }

        sourceField /= this->db().time().deltaTValue() * this->mesh().V();
    }

    return trhoTrans;
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingCloudSTFS<CloudType>::Srho(volScalarField& rho) const
{
    if (this->solution().coupled())
    {
        tmp<volScalarField> trhoTrans(
            new volScalarField(
                IOobject(
                    this->name() + ":rhoTrans",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false),
                this->mesh(),
                dimensionedScalar(dimMass / dimTime / dimVolume, Zero)));

        scalarField& sourceField = trhoTrans.ref();

        if (this->solution().semiImplicit("rho"))
        {

            forAll(rhoTrans_, i)
            {
                sourceField += rhoTrans_[i];
            }
            sourceField /= this->db().time().deltaTValue() * this->mesh().V();

            return fvm::SuSp(trhoTrans() / rho, rho);
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass / dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            forAll(rhoTrans_, i)
            {
                sourceField += rhoTrans_[i];
            }

            fvm.source() = -trhoTrans() / this->db().time().deltaT();

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>::New(rho, dimMass / dimTime);
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloudSTFS<CloudType>::haTrans()
{
    return *haTrans_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloudSTFS<CloudType>::haTrans() const
{
    return *haTrans_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloudSTFS<CloudType>::haCoeff()
{
    return *haCoeff_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloudSTFS<CloudType>::haCoeff() const
{
    return *haCoeff_;
}

template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingCloudSTFS<CloudType>::Sha(volScalarField& ha) const
{
    // if (debug)
    //{
    Info << "haTrans min/max = " << min(haTrans()).value() << ", "
         << max(haTrans()).value() << nl
         << "haCoeff min/max = " << min(haCoeff()).value() << ", "
         << max(haCoeff()).value() << endl;
    //}

    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("h"))
        {
            const volScalarField Cp(this->thermo().thermo().Cp());
            const volScalarField::Internal
                Vdt(this->mesh().V() * this->db().time().deltaT());

            return haTrans() / Vdt
                 - fvm::SuSp(haCoeff() / (Cp * Vdt), ha)
                 + haCoeff() / (Cp * Vdt) * ha;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(ha, dimEnergy / dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            fvm.source() = -haTrans() / (this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>::New(ha, dimEnergy / dimTime);
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloudSTFS<CloudType>::KTrans()
{
    return *KTrans_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloudSTFS<CloudType>::KTrans() const
{
    return *KTrans_;
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingCloudSTFS<CloudType>::SK(volScalarField& K) const
{
    if (debug)
    {
        Info << "KTrans min/max = " << min(KTrans()).value() << ", "
             << max(KTrans()).value() << endl;
    }

    if (this->solution().coupled())
    {
        // if (this->solution().semiImplicit("h"))
        //{
        //     const volScalarField Cp(thermo_.thermo().Cp());
        //     const volScalarField::Internal
        //         Vdt(this->mesh().V()*this->db().time().deltaT());

        //    return
        //        haTrans()/Vdt
        //      - fvm::SuSp(haCoeff()/(Cp*Vdt), ha)
        //      + haCoeff()/(Cp*Vdt)*ha;
        //}
        // else
        //{
        tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(K, dimEnergy / dimTime));
        fvScalarMatrix& fvm = tfvm.ref();

        fvm.source() = -KTrans() / (this->db().time().deltaT());

        return tfvm;
        //}
    }

    return tmp<fvScalarMatrix>::New(K, dimEnergy / dimTime);
}

template<class CloudType>
inline Foam::scalar
Foam::ReactingCloudSTFS<CloudType>::absoluteEnthalpyOfSystem() const
{
    scalar absoluteEnthalpy = 0;

    for (const parcelType& p : *this)
    {
        absoluteEnthalpy += p.nParticle() * p.mass() * p.ha();

        if (p.nParticle() < 0)
        {
            Info << "nParticle is negative: " << p.nParticle() << endl;
        }
        if (p.mass() < 0)
        {
            Info << "particle mass is negative: " << p.mass() << endl;
        }
        if (abs(p.ha() > 1e9))
        {
            Info << "absolute enthalpy of particle too great " << p.ha() << endl;
            Info << "ha particle too great, p.T, p.T, p.d  " << p.T() << " " << p.T() << " " << p.d() << endl;
        }
    }


    return absoluteEnthalpy;
}
template<class CloudType>
inline Foam::scalar Foam::ReactingCloudSTFS<CloudType>::Dmin() const
{
    scalar d = GREAT;
    for (const parcelType& p : *this)
    {
        d = min(d, p.d());
    }
    reduce(d, minOp<scalar>());
    return d;
}


// ************************************************************************* //
