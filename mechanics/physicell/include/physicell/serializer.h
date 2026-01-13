#pragma once

#include <memory>

#include <common/types.h>

namespace mechanics::physicell {

class serializer
{
public:
	virtual void serialize(real_t current_time) = 0;

	virtual ~serializer() = default;
};

using serializer_ptr = std::unique_ptr<serializer>;

} // namespace mechanics::physicell
