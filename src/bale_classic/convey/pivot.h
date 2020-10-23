// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// This file is designed to be "encased" by cases.pl.  The token Item
// must be #defcased to a number: either a positive item size or zero.

#if Item == 0
#define porter_push_1_0 porter_push
#define porter_push_2_0 porter_push
#define porter_push_4_0 porter_push
#endif

#defcases (Tag,Type) = (1,uint8_t), (4,uint32_t)
static bool
pivot_mid_##Tag##_##Item(tensor_t* matrix, buffer_t* buffer)
{
  // tag is (y'); we are (x',y); source is x, tag becomes (x)
  ACT_START(matrix_pivot);
  const uint64_t source = buffer->source;
  const size_t packet_bytes = (Item == 0) ? porter_stride(matrix->porters[0])
#if MATRIX_REMOTE_HOP == 0
    : Tag * (1 + (Item + Tag - 1) / Tag);
#else
    : 4 * (1 + (Item + 3) / 4);
#endif

  char* packet = buffer->data + buffer->start;
  char* limit = buffer->data + buffer->limit;
  for (; packet < limit; packet += packet_bytes) {
#if MATRIX_REMOTE_HOP == 0
    int dest = *(Type*)packet;
    bool ok = porter_push_4_##Item
      (matrix->porters[1], source, packet + Tag, dest);
#else
    int dest = *(uint32_t*)packet;
    bool ok = porter_push_##Tag##_##Item
      (matrix->porters[1], source, packet + 4, dest);
#endif
    if (!ok) {
      buffer->start = packet - buffer->data;
      ACT_STOP(matrix_pivot);
      return false;
    }
  }

  buffer->start = buffer->limit;
  ACT_STOP(matrix_pivot);
  return true;
}
#endcases

#defcases Tag = 1, 2, 4
static bool
pivot_early_##Tag##_##Item(tensor_t* tensor, buffer_t* buffer)
{
  // tag is (x',z'); source is z; we are (x,y,y')
  // hop to (x',y,y'), tag becomes (z, z') [4 bits + 4 bits if Tag==1]
  ACT_START(tensor_early);
  const uint32_t source = buffer->source << ((Tag == 1) ? 4 : 8);
  const size_t packet_bytes = (Item == 0) ?
    porter_stride(tensor->porters[0]) : 4 * (1 + (Item + 3) / 4);
  char* packet = buffer->data + buffer->start;
  char* limit = buffer->data + buffer->limit;
  for (; packet < limit; packet += packet_bytes) {
    uint32_t tag = *(uint32_t*)packet;
    int dest = tag >> 8;
    tag = (tag & 0xFF) | source;
    bool ok = porter_push_##Tag##_##Item
      (tensor->porters[1], tag, packet + 4, dest);
    if (!ok) {
      buffer->start = packet - buffer->data;
      ACT_STOP(tensor_early);
      return false;
    }
  }
  buffer->start = buffer->limit;
  ACT_STOP(tensor_early);
  return true;
}
#endcases

#defcases (Tag,Type) = (1,uint8_t), (2,uint16_t), (4,uint32_t)
static bool
pivot_late_##Tag##_##Item(tensor_t* tensor, buffer_t* buffer)
{
  // tag is (z,z'), source is x, we are (x',y',y)
  // hop to (x',y',z'), tag becomes (x,z)
  ACT_START(tensor_late);
  const uint32_t source = buffer->source << 8;
  const size_t packet_bytes = (Item == 0) ?
    porter_stride(tensor->porters[1]) : Tag * (1 + (Item + Tag - 1) / Tag);
  char* packet = buffer->data + buffer->start;
  char* limit = buffer->data + buffer->limit;
  for (; packet < limit; packet += packet_bytes) {
    uint32_t tag = *(Type*)packet;
    int dest = tag & ((Tag == 1) ? 0xF : 0xFF);
    tag = (tag >> ((Tag == 1) ? 4 : 8)) | source;
    bool ok = porter_push_4_##Item
      (tensor->porters[2], tag, packet + Tag, dest);
    if (!ok) {
      buffer->start = packet - buffer->data;
      ACT_STOP(tensor_late);
      return false;
    }
  }
  buffer->start = buffer->limit;
  ACT_STOP(tensor_late);
  return true;
}
#endcases
