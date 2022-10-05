#include "vertex_appearance.h"
bool operator==(const VertexAppearance& lhs, const VertexAppearance& rhs)
{
  return (lhs.v == rhs.v) && (lhs.time == rhs.time);
}
