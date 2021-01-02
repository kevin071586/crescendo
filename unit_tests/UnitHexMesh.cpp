#include <gtest/gtest.h>
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"

TEST(HexMesh, ElemId) {
    // Test 3x1x1 HexFixture structure
    const unsigned NX = 3;
    const unsigned NY = 1;
    const unsigned NZ = 1;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
    hf.m_meta.commit();
    hf.generate_mesh();
    EXPECT_EQ( hf.elem_id(0,0,0), 1u );
    EXPECT_EQ( hf.elem_id(1,0,0), 2u );
    EXPECT_EQ( hf.elem_id(2,0,0), 3u );
}

TEST(HexMesh, HasCoordinateField) {
    // Test 3x1x1 HexFixture structure
    const unsigned NX = 3;
    const unsigned NY = 1;
    const unsigned NZ = 1;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
    hf.m_meta.commit();
    hf.generate_mesh();

    auto coordField = hf.m_meta.get_field(stk::topology::NODE_RANK, "coordinates");
    EXPECT_FALSE(coordField == nullptr);
}
