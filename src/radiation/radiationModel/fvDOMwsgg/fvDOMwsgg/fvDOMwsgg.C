/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "fvDOMwsgg.H"
#include "absorptionEmissionModel.H"
#include "wsggAbsorptionEmission.H"
#include "wsggBinaryAbsorptionEmission.H"
#include "scatterModel.H"
#include "constants.H"
#include "unitConversion.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// #include "STFScoalCloud.H"


using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
defineTypeNameAndDebug(fvDOMwsgg, 0);
addToRadiationRunTimeSelectionTables(fvDOMwsgg);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::fvDOMwsgg::rotateInitialRays(const vector& sunDir)
{
    // Rotate Y spherical cordinates to Sun direction.
    // Solid angles on the equator are better fit for planar radiation
    const tensor coordRot = rotationTensor(vector(0, 1, 0), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }
}


void Foam::radiation::fvDOMwsgg::alignClosestRayToSun(const vector& sunDir)
{
    label SunRayId(-1);
    scalar maxSunRay = -GREAT;

    // Looking for the ray closest to the Sun direction
    forAll(IRay_, rayId)
    {
        const vector& iD = IRay_[rayId].d();
        scalar dir = sunDir & iD;
        if (dir > maxSunRay)
        {
            maxSunRay = dir;
            SunRayId = rayId;
        }
    }

    // Second rotation to align colimated radiation with the closest ray
    const tensor coordRot = rotationTensor(IRay_[SunRayId].d(), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }

    Info << "Sun direction : " << sunDir << nl << endl;
    Info << "Sun ray ID : " << SunRayId << nl << endl;
}


void Foam::radiation::fvDOMwsgg::updateRaysDir()
{
    solarCalculator_->correctSunDirection();
    const vector sunDir = solarCalculator_->direction();

    // First iteration
    if (updateTimeIndex_ == 0)
    {
        rotateInitialRays(sunDir);
        alignClosestRayToSun(sunDir);
    }
    else if (updateTimeIndex_ > 0)
    {
        alignClosestRayToSun(sunDir);
    }
}


void Foam::radiation::fvDOMwsgg::initialise()
{
    coeffs_.readIfPresent("useExternalBeam", useExternalBeam_);

    if (useExternalBeam_)
    {
        coeffs_.readEntry("spectralDistribution", spectralDistribution_);

        spectralDistribution_ =
            spectralDistribution_ / sum(spectralDistribution_);

        const dictionary& solarDict = this->subDict("solarCalculatorCoeffs");
        solarCalculator_.reset(new solarCalculator(solarDict, mesh_));

        if (mesh_.nSolutionD() != 3)
        {
            FatalErrorInFunction
                << "External beam model only available in 3D meshes "
                << abort(FatalError);
        }

        if (solarCalculator_->diffuseSolarRad() > 0)
        {
            FatalErrorInFunction
                << "External beam model does not support Diffuse "
                << "Solar Radiation. Set diffuseSolarRad to zero"
                << abort(FatalError);
        }
        if (spectralDistribution_.size() != nLambda_)
        {
            FatalErrorInFunction
                << "The epectral energy distribution has different bands "
                << "than the absoprtivity model "
                << abort(FatalError);
        }
    }

    // 3D
    if (mesh_.nSolutionD() == 3)
    {
        nRay_ = 4 * nPhi_ * nTheta_;

        IRay_.setSize(nRay_);

        const scalar deltaPhi = pi / (2 * nPhi_);
        const scalar deltaTheta = pi / nTheta_;

        label i = 0;

        for (label n = 1; n <= nTheta_; n++)
        {
            for (label m = 1; m <= 4 * nPhi_; m++)
            {
                scalar thetai = (2 * n - 1) * deltaTheta / 2.0;
                scalar phii = (2 * m - 1) * deltaPhi / 2.0;

                IRay_.set(
                    i,
                    new radiativeIntensityRaywsgg(
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        *absorptionEmission_,
                        blackBody_,
                        i));
                i++;
            }
        }
    }
    // 2D
    else if (mesh_.nSolutionD() == 2)
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;
        nRay_ = 4 * nPhi_;
        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi / (2.0 * nPhi_);
        label i = 0;
        for (label m = 1; m <= 4 * nPhi_; m++)
        {
            const scalar phii = (2 * m - 1) * deltaPhi / 2.0;
            IRay_.set(
                i,
                new radiativeIntensityRaywsgg(
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    *absorptionEmission_,
                    blackBody_,
                    i));
            i++;
        }
    }
    // 1D
    else
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;
        nRay_ = 2;
        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi;
        label i = 0;
        for (label m = 1; m <= 2; m++)
        {
            const scalar phii = (2 * m - 1) * deltaPhi / 2.0;
            IRay_.set(
                i,
                new radiativeIntensityRaywsgg(
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    *absorptionEmission_,
                    blackBody_,
                    i));
            i++;
        }
    }


    // Construct absorption field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set(
            lambdaI,
            new volScalarField(
                IOobject(
                    "aLambda_" + Foam::name(lambdaI),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                a_));
    }

    Info << "fvDOMwsgg : Allocated " << IRay_.size()
         << " rays with average orientation:" << nl;

    if (useExternalBeam_)
    {
        // Rotate rays for Sun direction
        updateRaysDir();
    }

    scalar totalOmega = 0;
    forAll(IRay_, rayId)
    {
        if (omegaMax_ < IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        totalOmega += IRay_[rayId].omega();
        Info << '\t' << IRay_[rayId].I().name() << " : " << "dAve : "
             << '\t' << IRay_[rayId].dAve() << " : " << "omega : "
             << '\t' << IRay_[rayId].omega() << " : " << "d : "
             << '\t' << IRay_[rayId].d() << nl;
    }

    Info << "Total omega : " << totalOmega << endl;

    Info << endl;

    // Give the thermo pointer to the absorption/emission model, if it needs one
    // (wsggAbsorptionEmission-derived models resolve their thermo lazily, since
    // it isn't available at their own construction time).
    if (isA<wsggAbsorptionEmission>(absorptionEmission_()))
    {
        dynamic_cast<wsggAbsorptionEmission&>(absorptionEmission_())
            .setThermoPtr(*thermoPtr());
    }
    else if (isA<wsggBinaryAbsorptionEmission>(absorptionEmission_()))
    {
        dynamic_cast<wsggBinaryAbsorptionEmission&>(absorptionEmission_())
            .setThermoPtr(*thermoPtr());
        Info << "detected wsggBinaryAbsorptionEmission model" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fvDOMwsgg::fvDOMwsgg(const volScalarField& T)
: radiationModelSTFS(typeName, T),
  G_(
      IOobject(
          "G",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar(dimMass / pow3(dimTime), Zero)),
  //
  kG_(
      IOobject(
          "kG",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar("kG", dimMass / pow3(dimTime) / dimLength, 0.0)),
  //
  qr_(
      IOobject(
          "qr",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar(dimMass / pow3(dimTime), Zero)),
  qem_(
      IOobject(
          "qem",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE),
      mesh_,
      dimensionedScalar(dimMass / pow3(dimTime), Zero)),
  qin_(
      IOobject(
          "qin",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar(dimMass / pow3(dimTime), Zero)),
  a_(
      IOobject(
          "a",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE),
      mesh_,
      dimensionedScalar(dimless / dimLength, Zero)),

  nTheta_(coeffs_.get<label>("nTheta")),
  nPhi_(coeffs_.get<label>("nPhi")),
  nRay_(0),
  nLambda_(absorptionEmission_->nBands()),
  aLambda_(nLambda_),
  blackBody_(nLambda_, T),
  IRay_(0),
  tolerance_(
      coeffs_.getOrDefaultCompat<scalar>(
          "tolerance",
          {{"convergence", 1712}},
          0)),
  maxIter_(coeffs_.getOrDefault<label>("maxIter", 50)),
  omegaMax_(0),
  meshOrientation_(
      coeffs_.getOrDefault<vector>("meshOrientation", Zero)),
  useExternalBeam_(false),
  spectralDistribution_(),
  solarCalculator_(),
  updateTimeIndex_(0)
{
    initialise();
}


Foam::radiation::fvDOMwsgg::fvDOMwsgg(
    const dictionary& dict,
    const volScalarField& T)
: radiationModelSTFS(typeName, dict, T),
  G_(
      IOobject(
          "G",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar(dimMass / pow3(dimTime), Zero)),
  //
  kG_(
      IOobject(
          "kG",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar("kG", dimMass / pow3(dimTime) / dimLength, 0.0)),
  //
  qr_(
      IOobject(
          "qr",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar(dimMass / pow3(dimTime), Zero)),
  qem_(
      IOobject(
          "qem",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE),
      mesh_,
      dimensionedScalar(dimMass / pow3(dimTime), Zero)),
  qin_(
      IOobject(
          "qin",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar(dimMass / pow3(dimTime), Zero)),
  a_(
      IOobject(
          "a",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE),
      mesh_,
      dimensionedScalar(dimless / dimLength, Zero)),

  nTheta_(coeffs_.get<label>("nTheta")),
  nPhi_(coeffs_.get<label>("nPhi")),
  nRay_(0),
  nLambda_(absorptionEmission_->nBands()),
  aLambda_(nLambda_),
  blackBody_(nLambda_, T),
  IRay_(0),
  tolerance_(
      coeffs_.getOrDefaultCompat<scalar>(
          "tolerance",
          {{"convergence", 1712}},
          0)),
  maxIter_(coeffs_.getOrDefault<label>("maxIter", 50)),
  omegaMax_(0),
  meshOrientation_(
      coeffs_.getOrDefault<vector>("meshOrientation", Zero)),
  useExternalBeam_(false),
  spectralDistribution_(),
  solarCalculator_(),
  updateTimeIndex_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::fvDOMwsgg::~fvDOMwsgg()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::fvDOMwsgg::read()
{
    if (radiationModelSTFS::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresentCompat(
            "tolerance", {{"convergence", 1712}}, tolerance_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }

    return false;
}


void Foam::radiation::fvDOMwsgg::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);
    // Foam::tmp<Foam::volScalarField> e = absorptionEmission_->eCont();
    updateBlackBodyEmission();

    if (useExternalBeam_)
    {
        switch (solarCalculator_->sunDirectionModel())
        {
        case solarCalculator::mSunDirConstant:
        {
            break;
        }
        case solarCalculator::mSunDirTracking:
        {
            label updateIndex = label(
                mesh_.time().value()
                / solarCalculator_->sunTrackingUpdateInterval());

            if (updateIndex > updateTimeIndex_)
            {
                Info << "Updating Sun position..." << endl;
                updateTimeIndex_ = updateIndex;
                updateRaysDir();
            }
            break;
        }
        }
    }

    // Set rays convergence false
    List<bool> rayIdConv(nRay_, false);

    scalar maxResidual = 0;
    label radIter = 0;
    do
    {
        Info << "Radiation solver iter: " << radIter << endl;

        radIter++;
        maxResidual = 0;
        forAll(IRay_, rayI)
        {
            if (!rayIdConv[rayI])
            {
                scalar maxBandResidual = IRay_[rayI].correct();
                maxResidual = max(maxBandResidual, maxResidual);

                if (maxBandResidual < tolerance_)
                {
                    rayIdConv[rayI] = true;
                }
            }
        }

    } while (maxResidual > tolerance_ && radIter < maxIter_);

    updateG();

    // Here a write function can be added for source term
}


Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOMwsgg::Rp() const
{
    // Construct using contribution from first frequency band
    tmp<volScalarField> tRp(
        new volScalarField(
            IOobject(
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false),
            (
                4
                * physicoChemical::sigma
                * (absorptionEmission_->eCont(0)) * (absorptionEmission_->aCont(0)))));

    volScalarField& Rp = tRp.ref();

    // Add contributions over remaining frequency bands
    for (label j = 1; j < nLambda_; j++)
    {
        Rp +=
            (4
             * physicoChemical::sigma
             * absorptionEmission_->eCont(j) * absorptionEmission_->aCont(j));
    }

    if (db().time().outputTime())
    {
        Rp.write();
    }

    return tRp;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::fvDOMwsgg::Ru() const
{
    tmp<DimensionedField<scalar, volMesh>> tRu(
        new DimensionedField<scalar, volMesh>(
            IOobject(
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false),
            mesh_,
            dimensionedScalar(dimensionSet(1, -1, -3, 0, 0), Zero)));

    DimensionedField<scalar, volMesh> Eprt =
        absorptionEmission_->EDisp(0);

    for (label j = 1; j < nLambda_; j++)
    {
        Eprt += absorptionEmission_->EDisp(j);
    }


    DimensionedField<scalar, volMesh>& Ru = tRu.ref();

    const DimensionedField<scalar, volMesh> E =
        absorptionEmission_->ECont()()();

    Ru = kG_; // - E;
              //	-Eprt;  //Ru and Rp are only for gas phase, particle source is calculated in lagrangian class

    ////////////////////////////////////////////////
    // Debug particles!

    //    const kinematicCloud& kcloud =
    //        mesh_.objectRegistry::lookupObject<kinematicCloud>("coalCloud"); //Name should be checked properly
    //
    //    const STFScoalCloud& cloud =
    //        dynamic_cast<const STFScoalCloud&>(kcloud);
    //
    //          //tmp<DimensionedField<scalar, volMesh>> Shradp
    //          //(new DimensionedField<scalar, volMesh>
    //          volScalarField Shradp
    //	  //(volScalarField
    //            (
    //                IOobject
    //                (
    //                    "Shradp",
    //                    mesh_.time().timeName(),
    //                    mesh_,
    //                    IOobject::NO_READ,
    //                    IOobject::AUTO_WRITE
    //                ),
    //                mesh_,
    //                dimensionedScalar("Shradp", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    //                //zeroGradientFvPatchScalarField::typeName
    //            )	;
    //         // );
    //
    //        DimensionedField<scalar, volMesh>&  ShradPrt = Shradp.ref();
    //
    //	scalar cloudEm_ = cloud.constProps().epsilon0();
    //    for(const auto& c : cloud)
    //    {
    //	ShradPrt[c.cell()] += //c.m()*c.Cp()*c.deltaT/dt
    //        		    c.areaS() * cloudEm_ * c.nParticle() / mesh_.V()[c.cell()]
    //			    * 	(
    //			     	G_[c.cell()]/4.0
    //				 - constant::physicoChemical::sigma.value()*pow4(c.T())
    //				);
    //			    //	);
    ////			     c.areaS() * cloudEm_ * c.nParticle() / mesh_.V()[c.cell()]
    ////				*( - constant::physicoChemical::sigma.value()*pow4(c.T()));
    //				//-4.0*c.Eb();
    //    }
    //    //ShradPrt = ShradPrt + Eprt;
    //    if(db().time().outputTime())
    //    {
    //	ShradPrt.write();
    //    }
    //
    //
    //
    //////////////////////////////////////

    if (db().time().outputTime())
    {
        Ru.write();
    }
    return tRu;
}


void Foam::radiation::fvDOMwsgg::updateBlackBodyEmission()
{
    for (label j = 0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}


void Foam::radiation::fvDOMwsgg::updateG()
{
    G_ = dimensionedScalar(dimMass / pow3(dimTime), Zero);
    kG_ = dimensionedScalar("zero", dimMass / pow3(dimTime) / dimLength, 0.0);
    qr_ = dimensionedScalar(dimMass / pow3(dimTime), Zero);
    qem_ = dimensionedScalar(dimMass / pow3(dimTime), Zero);
    qin_ = dimensionedScalar(dimMass / pow3(dimTime), Zero);

    forAll(IRay_, rayI)
    {
        IRay_[rayI].addIntensity();
        volScalarField temp = IRay_[rayI].I();
        G_ += IRay_[rayI].I() * IRay_[rayI].omega();
        kG_ += IRay_[rayI].kI() * IRay_[rayI].omega();

        qr_.boundaryFieldRef() += IRay_[rayI].qr().boundaryField();
        qem_.boundaryFieldRef() += IRay_[rayI].qem().boundaryField();
        qin_.boundaryFieldRef() += IRay_[rayI].qin().boundaryField();
    }

    // kG_.write();

    // qem_.write();
    // qin_.write();
}


void Foam::radiation::fvDOMwsgg::setRayIdLambdaId(
    const word& name,
    label& rayId,
    label& lambdaId) const
{
    // Assuming name is in the form: CHARS_rayId_lambdaId
    const auto i1 = name.find('_');
    const auto i2 = name.find('_', i1 + 1);

    rayId = readLabel(name.substr(i1 + 1, i2 - i1 - 1));
    lambdaId = readLabel(name.substr(i2 + 1));
}


const Foam::solarCalculator& Foam::radiation::fvDOMwsgg::solarCalc() const
{
    return solarCalculator_();
}


// ************************************************************************* //
