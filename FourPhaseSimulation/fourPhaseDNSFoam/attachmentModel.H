/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/ and interFoam/interMixingFoam/ (solver family)
    Copyright (C) OpenFOAM Foundation and OpenCFD Ltd.; see the upstream
    file for the exact original copyright years.
    Developed and ported mainly against OpenFOAM v2512; the exact upstream
    version originally copied has not been verified.

Modifications
    Copyright (C) 2026 Bjoern Klepsch
    Last modified: 2026-07-20 (see git history for the detailed change log)
-------------------------------------------------------------------------------
License
    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See the LICENSE file in the repository root.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
\*---------------------------------------------------------------------------*/

#ifndef attachmentModel_H
#define attachmentModel_H

#include "fvCFD.H"
#include "basicKinematicCloud.H"
#include "mathematicalConstants.H"
#include "Random.H"

class attachmentModel
{
    scalar Pc_;
    scalar Pa_;
    scalar kDetach_;
    scalar rCapture_;

    List<Tuple2<label,label>> attachedPairs_;
    mutable Random rng_;

public:

    attachmentModel(const dictionary& dict)
    :
        Pc_(dict.get<scalar>("collisionProbability")),
        Pa_(dict.get<scalar>("attachmentEfficiency")),
        kDetach_(dict.get<scalar>("detachmentRate")),
        rCapture_(dict.get<scalar>("captureRadius")),
        rng_(Foam::clock::getTime())
    {}

    void solve
    (
        basicKinematicCloud& particles,
        basicKinematicCloud& bubbles,
        scalar dt
    )
    {
        forAllIter(basicKinematicCloud, particles, pIter)
        {
            auto& p = pIter();
            if (p.active())
            {
                forAllIter(basicKinematicCloud, bubbles, bIter)
                {
                    auto& b = bIter();
                    if (b.active())
                    {
                        scalar dist = mag(p.position() - b.position());
                        if (dist < rCapture_)
                        {
                            scalar P_attach = Pc_ * Pa_;
                            scalar randNum = rng_.sample01<scalar>();

                            if (randNum < P_attach)
                            {
                                scalar mp = p.mass();
                                scalar mb = b.mass();
                                scalar Vp = mp / p.rho();
                                scalar Vb = mb / b.rho();
                                scalar rhoAgg = (mp + mb) / (Vp + Vb);

                                b.U() = (mp*p.U() + mb*b.U()) / (mp + mb);
                                // p.active() = false  // NOTE: use cloud remove API instead;
                                p.U() = b.U();

                                attachedPairs_.append
                                (
                                    Tuple2<label,label>(p.origId(), b.origId())
                                );

                                Info<< "Attachment: particle " << p.origId()
                                    << " attached to bubble " << b.origId()
                                    << " | rho_agg = " << rhoAgg << endl;
                            }
                        }
                    }
                }
            }
        }

        // Detachment
        scalar Pdet = 1.0 - Foam::exp(-kDetach_ * dt);

        for (label i = attachedPairs_.size()-1; i >= 0; i--)
        {
            scalar r = rng_.sample01<scalar>();
            if (r < Pdet)
            {
                label pId = attachedPairs_[i].first();
                label bId = attachedPairs_[i].second();

                forAllIter(basicKinematicCloud, particles, pIter)
                {
                    if (pIter().origId() == pId)
                    {
                        // pIter().active() = true; // read-only in OF2512
                        Info<< "Detachment: particle " << pId
                            << " released from bubble " << bId << endl;
                    }
                }
                attachedPairs_[i] = attachedPairs_.last(); attachedPairs_.resize(attachedPairs_.size()-1);
            }
        }
    }
};

#endif
