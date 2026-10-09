// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. See the COPYRIGHT file at the top level of the distribution.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_IBAngularSpringForceSpec
#define included_IBAMR_IBAngularSpringForceSpec

#include <ibamr/config.h>

#include <ibtk/Streamable.h>
#include <ibtk/StreamableFactory.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <utility>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class IntVector;
} // namespace hier
namespace tbox
{
class AbstractStream;
} // namespace tbox
} // namespace SAMRAI

namespace IBAMR
{
/*! \brief Data for angular springs attached to one Lagrangian node.
 *
 * An angular spring connects a master node to two neighbors. Its energy is
 * \f$\tfrac{1}{2} k\theta^2\f$, where \f$\theta\f$ is the turning angle.
 */
class IBAngularSpringForceSpec : public IBTK::Streamable
{
public:
    using NeighborIdxs = std::pair<int, int>;

    static void registerWithStreamableManager();
    static bool getIsRegisteredWithStreamableManager();
    static int STREAMABLE_CLASS_ID;

    IBAngularSpringForceSpec(unsigned int num_springs = 0);
    IBAngularSpringForceSpec(int master_idx,
                             const std::vector<NeighborIdxs>& neighbor_idxs,
                             const std::vector<double>& stiffnesses,
                             const std::vector<double>& rest_angles);
    ~IBAngularSpringForceSpec();

    unsigned int getNumberOfAngularSprings() const;
    const int& getMasterNodeIndex() const;
    int& getMasterNodeIndex();
    const std::vector<NeighborIdxs>& getNeighborNodeIndices() const;
    std::vector<NeighborIdxs>& getNeighborNodeIndices();
    const std::vector<double>& getStiffnesses() const;
    std::vector<double>& getStiffnesses();
    const std::vector<double>& getRestAngles() const;
    std::vector<double>& getRestAngles();
    int getStreamableClassID() const override;
    size_t getDataStreamSize() const override;
    void packStream(SAMRAI::tbox::AbstractStream& stream) override;

private:
    IBAngularSpringForceSpec(const IBAngularSpringForceSpec&) = delete;
    IBAngularSpringForceSpec& operator=(const IBAngularSpringForceSpec&) = delete;

    int d_master_idx = IBTK::invalid_index;
    std::vector<NeighborIdxs> d_neighbor_idxs;
    std::vector<double> d_stiffnesses;
    std::vector<double> d_rest_angles;

    class Factory : public IBTK::StreamableFactory
    {
    public:
        ~Factory() = default;
        int getStreamableClassID() const override;
        void setStreamableClassID(int class_id) override;
        SAMRAI::tbox::Pointer<IBTK::Streamable> unpackStream(SAMRAI::tbox::AbstractStream& stream,
                                                             const SAMRAI::hier::IntVector<NDIM>& offset) override;

    private:
        Factory();
        Factory(const Factory&) = delete;
        Factory& operator=(const Factory&) = delete;
        friend class IBAngularSpringForceSpec;
    };
    using IBAngularSpringForceSpecFactory = IBAngularSpringForceSpec::Factory;
};
} // namespace IBAMR

#include <ibamr/private/IBAngularSpringForceSpec-inl.h> // IWYU pragma: keep

#endif // #ifndef included_IBAMR_IBAngularSpringForceSpec
