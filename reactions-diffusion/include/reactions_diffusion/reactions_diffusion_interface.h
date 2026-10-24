#pragma once

#include <array>
#include <memory>
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

	virtual real_t get_substrate_density(index_t s, std::span<const real_t> position) const = 0;

	virtual std::array<real_t, 3> get_substrate_gradient(index_t s, std::span<const real_t> position) const = 0;

	~reactions_diffusion_interface() override = default;
};

using reactions_diffusion_interface_ptr = std::shared_ptr<reactions_diffusion_interface>;

} // namespace physicore::reactions_diffusion
