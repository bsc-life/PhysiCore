#pragma once

#include <memory>
#include <optional>

#include <common/timestep_executor.h>
#include <common/types.h>
#include <common/cartesian_mesh.h>

#include "mechanical_agent_container.h"
#include "serializer.h"
#include "solver.h"

namespace physicore::mechanics::physicell {

class environment : public timestep_executor
{
public:
	environment(real_t timestep, index_t dims, index_t agent_types_count, index_t substrates_count);

	void run_single_timestep() override;

	void serialize_state(real_t current_time) override;


	real_t timestep;
	bool automated_spring_adhesion;
	bool virtual_wall_at_domain_edges;

	serializer_ptr serializer;
	solver_ptr solver;

	std::unique_ptr<mechanical_agent_container> agents;

	// Mesh data for spatial queries 
	void set_mesh(cartesian_mesh mesh) { this->mesh_ = std::move(mesh); }
	const cartesian_mesh& get_mesh() const;
	bool has_mesh() const { return mesh_.has_value();}
private:
	std::optional<cartesian_mesh> mesh_;
};

} // namespace physicore::mechanics::physicell
