#pragma once

#include <span>
#include <string>

#include <common/timestep_executor.h>
#include <common/types.h>

namespace physicore::reactions_diffusion {

class reactions_diffusion_interface : public timestep_executor
{
public:
	reactions_diffusion_interface() = default;
	reactions_diffusion_interface(const reactions_diffusion_interface&) = delete;
	reactions_diffusion_interface(reactions_diffusion_interface&&) = delete;
	reactions_diffusion_interface& operator=(const reactions_diffusion_interface&) = delete;
	reactions_diffusion_interface& operator=(reactions_diffusion_interface&&) = delete;

	virtual std::span<const std::string> get_substrate_names() const = 0;
	virtual std::span<const std::string> get_substrate_units() const = 0;

	virtual real_t get_substrate_density(index_t s, index_t x, index_t y, index_t z) const = 0;

	~reactions_diffusion_interface() override = default;
};

} // namespace physicore::reactions_diffusion
