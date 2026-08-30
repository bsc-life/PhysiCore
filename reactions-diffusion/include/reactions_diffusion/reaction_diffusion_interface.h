#pragma once

#include <span>
#include <string>

#include <common/timestep_executor.h>
#include <common/types.h>

namespace physicore::reactions_diffusion {

class reaction_diffusion_interface : public timestep_executor
{
public:
	reaction_diffusion_interface() = default;
	reaction_diffusion_interface(const reaction_diffusion_interface&) = delete;
	reaction_diffusion_interface(reaction_diffusion_interface&&) = delete;
	reaction_diffusion_interface& operator=(const reaction_diffusion_interface&) = delete;
	reaction_diffusion_interface& operator=(reaction_diffusion_interface&&) = delete;

	virtual std::span<const std::string> get_substrate_names() const = 0;
	virtual std::span<const std::string> get_substrate_units() const = 0;

	virtual real_t get_substrate_density(index_t s, index_t x, index_t y, index_t z) const = 0;

	~reaction_diffusion_interface() override = default;
};

} // namespace physicore::reactions_diffusion
